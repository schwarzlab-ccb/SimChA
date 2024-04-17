using System.Diagnostics;
using System.Globalization;
using CommandLine;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Optimization;
using SimChA.Simulation;

Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
var cmdOptions = Parser.Default.ParseArguments<CmdOptions>(args);
cmdOptions.WithNotParsed(o =>
{
    Console.WriteLine("Exiting");
    o.ToList().ForEach(Console.WriteLine); // Write out errors
    Environment.Exit(1);
});
var options = cmdOptions.Value;
var execMode = options.ExecMode;

SimParams simParams = FileIO.ReadSimParams(options.ConfigFile);

var rnd = new Random(simParams.Seed);
var files = new FileIO(options.OutputPath);
bool parseGenContents = execMode == ExecMode.ParseGenContents;
var includeSexChromosomes = !options.AutosomesOnly;
var genRef = FileIO.GetGenRef(options.DataFolder, includeSexChromosomes, parseGenContents);

var watch = new Stopwatch();
watch.Start();
List<Sample> samples;
if (execMode == ExecMode.BinSamples)
{
    if (options.CNProfiles == "")
    {
        throw new Exception("Error: No CN profiles provided. Cannot bin samples.");
    }
    Console.WriteLine("Reading observed data:");
    var profiles = FileIO.ReadProfiles(genRef, options.CNProfiles);
    var observedSamples = Simulator.SamplesFromProfiles(profiles);
    var binner = new Binner(genRef);
    Console.WriteLine("Binning samples:");
    var binnedSamples = binner.GetBinnedCNProfiles(observedSamples);
    Console.WriteLine("Writing binned samples:");
    files.WriteCopyNumbers(binnedSamples);
    return 0;
}
else if (execMode == ExecMode.RunOptimization)
{
    if (simParams.OptimizationParams is null)
    {
        throw new Exception("Error: OptimizationParams not set. Cannot perform optimization. Please set OptimizationParams.");
    }
    
    if (simParams.OptimizationParams.Mode == "Events")
    {
        Console.WriteLine("Event Optimization Mode -------- ");
        var optimizer = new Optimizer(simParams, rnd, options.Repeats, genRef, includeSexChromosomes);
        if (options.TargetParams != "")
        {
            Console.WriteLine("Generating Simulated Data");
            SimParams targetParams = FileIO.ReadSimParams(options.TargetParams);
            optimizer.InitializeObservations(targetParams);
        }
        else if (options.CNProfiles != "" && options.EventCounts != "")
        {
            Console.WriteLine("Reading observed data:");
            var profiles = FileIO.ReadProfiles(genRef, options.CNProfiles);
            var eventCounts = FileIO.ReadEventCounts(options.EventCounts);
            var observedSamples = Simulator.SamplesFromProfiles(profiles);
            optimizer.InitializeObservations(observedSamples, eventCounts);
        }
        else 
        {
            throw new Exception("Error: No target parameters (synthetic data) or bootstrap file (observed data) provided. Cannot perform event optimization without a target data set.");
        }
        
        if (simParams.OptimizationParams.UseABC)
        {
            var dist = optimizer.GetABCDistance();
            Console.WriteLine($"ABC distance: {dist}");
            return 0;
        }
        else
        {
            var outParams = optimizer.Optimize(files);
            files.WriteSimParams(outParams);
        }
    }
    else if (simParams.OptimizationParams.Mode == "Fitness")
    {
        Console.WriteLine("Fitness Optimization Mode -------- ");
        var optimizer = new FitnessOptimizer(simParams, rnd, options.Repeats, genRef, includeSexChromosomes);
        if (options.TargetParams != "")
        {
            Console.WriteLine("Generating Simulated Data");
            SimParams targetParams = FileIO.ReadSimParams(options.TargetParams);
            optimizer.InitializeObservations(targetParams);
        }
        else if (options.BinnedSamples != "")
        {
            var cloneComponents = FileIO.ReadCloneComponents(options.BootstrapFile);
            //optimizer.InitializeObservations(cloneComponents);
        }
        else
        {
            throw new Exception("Error: No target parameters (synthetic data) or bootstrap file (observed data) provided. Cannot perform fitness optimization without a target data set.");
        }
        
        var outParams = optimizer.Optimize(files);
        files.WriteSimParams(outParams);
    }
    else
    {
        throw new Exception("Error: Optimization mode not recognized. Please set OptimizationParams.Mode to either Events or Fitness.");
    }
    watch.Stop();
    //Console.WriteLine($"Total time: {TimeSpan.FromMilliseconds(watch.ElapsedMilliseconds)}");
    Console.WriteLine("Optimization finished");
    return 0;
}
if (execMode == ExecMode.Profiles)
{
    Console.WriteLine("Reading profiles:");
    var profiles = FileIO.ReadProfiles(genRef, options.CNProfiles);
    samples = Simulator.SamplesFromProfiles(profiles);
}
else
{
    if (simParams.Signatures is null || simParams.Signatures.Count == 0)
    {
        throw new Exception("No signatures were provided.");
    }
    Validators.ValidateSignatures(simParams.Signatures);
    Console.WriteLine("Computing mutations:");
    var simulator = new Simulator(rnd, genRef);
    if (execMode == ExecMode.Tree)
    {
        var inClones = FileIO.ReadClones(options.CloneTreeFile, options.UseMCMC);
        var (cnEventPs, mixture) = Converters.PropagateSigs(simParams.Signatures);
        string sampleName = Path.GetFileNameWithoutExtension(options.CloneTreeFile);
        var treeSample = new Sample(sampleName, Sampling.GetBinarySex(rnd, simParams.Sex), inClones, cnEventPs, mixture);
        samples = new List<Sample> {treeSample};
    }
    else
    {
        if (simParams.MCTarget is null && options.UseMCMC)
        {
            throw new Exception("Error: MCTarget not set. Cannot perform MC sampling. Please set MCTarget in the config file.");
        }
        if (execMode == ExecMode.Bootstrap)
        {
            var clonesList = FileIO.ReadFitnesses(options.BootstrapFile, simParams.Fitness);
            samples = Converters.MakeSamples(rnd, options.Repeats, simParams.EventCount, simParams.EventDist, simParams.Signatures, simParams.Sex, clonesList);
        }
        else
        {
            samples = Converters.MakeSamples(rnd, options.Repeats, simParams.EventCount, simParams.EventDist, simParams.Signatures, simParams.Sex, simParams.MCTarget);
        }
    }

    foreach (var sample in samples)
    {
        if (options.UseMCMC)
        {
            if (simParams.MCParams == null)
            {
                throw new Exception("Error: MCParams not set. Cannot perform MC sampling. Please set MCParams.");
            }

            simulator = new MCSimulator(rnd, genRef, simParams.Fitness, simParams.MCParams);
        }
        simulator.SampleEvents(sample);
    }
}
files.WriteSimParams(simParams);

Console.WriteLine("");

// Fitness data
Console.WriteLine("Computing clone stats:");
foreach (var sample in samples)
{
    int counter = 1;
    int total = sample.Clones.Count;
    foreach (var clone in sample.Clones)
    {
        Console.Write($"\rSample {sample.SampleId}. Clone {counter++}/{total}.".PadRight(80));
        sample.Stats[clone.CloneId] = CNProfile.GetCloneStats(sample, clone, genRef, simParams.Fitness, sample.Kars);
    }
}

try
{
    Console.WriteLine("");
    files.WriteSamples(samples);
    files.WriteCopyNumbers(genRef, samples);
    if (options.CalcConsistentCNs)
    {
        files.WriteConsistentCNs(genRef, samples);
    }
    files.WriteClones(samples);
    files.WriteKaryotypes(samples);
    if (samples.Any(s => s.EventDescs.Any()))
    {
        files.WriteEvents(samples);
    }
    if (options.UseVariants)
    {
        files.WriteVCF(genRef, samples);
    }
    if (options.WriteFasta)
    {
        files.WriteFasta(genRef, samples);
    }
}
catch (Exception e)
{
    Console.WriteLine($"Failed with exception {e.Message}. Stack: {e.StackTrace}");
    return e.HResult;
}

watch.Stop();
Console.WriteLine($"Total time: {TimeSpan.FromMilliseconds(watch.ElapsedMilliseconds)}");

return 0;
