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
List<Sample> samples = new();
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
    
    var optimizer = new Optimizer(simParams, rnd, options.Repeats, genRef, includeSexChromosomes, files);
    if (simParams.OptimizationParams.Mode == "Events")
    {
        Console.WriteLine("Event Optimization Mode -------- ");
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
        var outParams = optimizer.Optimize();
        files.WriteSimParams(outParams);
    }
    else if (simParams.OptimizationParams.Mode == "Fitness")
    {
        Console.WriteLine("Fitness Optimization Mode -------- ");
        optimizer = new FitnessOptimizer(simParams, rnd, options.Repeats, genRef, includeSexChromosomes, files);
        if (options.TargetParams != "")
        {
            Console.WriteLine("Generating Simulated Data");
            SimParams targetParams = FileIO.ReadSimParams(options.TargetParams);
            optimizer.InitializeObservations(targetParams);
        }
        else if (options.EventCounts != "" && options.CNProfiles != "")
        {
            Console.WriteLine("Reading observed data:");
            var profiles = FileIO.ReadProfiles(genRef, options.CNProfiles);
            var observedSamples = Simulator.SamplesFromProfiles(profiles);
            var eventCounts = FileIO.ReadEventCounts(options.EventCounts);
            foreach (var sample in observedSamples)
            {
                int counter = 1;
                int total = sample.Clones.Count;
                foreach (var clone in sample.Clones)
                {
                    Console.Write($"\rSample {sample.SampleId}. Clone {counter++}/{total}.".PadRight(80));
                    sample.Stats[clone.CloneId] = CNProfile.GetCloneStats(sample, clone, genRef, simParams.Fitness, sample.Kars);
                }
            }
            optimizer.InitializeObservations(observedSamples, eventCounts);
        }
        else
        {
            throw new Exception("Error: No target parameters (synthetic data) or bootstrap file (observed data) provided. Cannot perform fitness optimization without a target data set.");
        }
        
        var outParams = optimizer.Optimize();
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
    if (options.FitnessLandscape)
    {
        Console.WriteLine("Computing fitness landscape:");
        FitnessLandscape.GenerateFitnessLandscape(genRef, simParams, samples, files);
    }
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
    else if (execMode == ExecMode.UseMCMC)
    {
        if (simParams.MCParams is null)
        {
            throw new Exception("Error: MCParams not set. Cannot perform MCMC sampling. Please set MCParams in the config file.");
        }

        if (options.CNProfiles != "" && options.EventCounts != "")
        {
            var profiles = FileIO.ReadProfiles(genRef, options.CNProfiles);
            var eventCounts = FileIO.ReadEventCounts(options.EventCounts);
            var fitnessList = simulator.FitnessListFromSamples(simParams, profiles, eventCounts);
            samples = Converters.MakeSamples(rnd, options.Repeats, simParams.EventCount, simParams.EventDist, simParams.Signatures, simParams.Sex, fitnessList);
        }
        else
        {
            if (simParams.MCTarget is null)
            {
                throw new Exception("Error: MCTarget not set. Cannot perform MC sampling. Please set MCTarget in the config file.");
            }
            samples = Converters.MakeSamples(rnd, options.Repeats, simParams.EventCount, simParams.EventDist, simParams.Signatures, simParams.Sex, simParams.MCTarget);
        }
        simulator = new MCSimulator(rnd, genRef, simParams.Fitness, simParams.MCParams);
    }
    else
    {
        samples = Converters.MakeSamples(rnd, options.Repeats, simParams.EventCount, simParams.EventDist, simParams.Signatures, simParams.Sex, simParams.MCTarget);
    }

    foreach (var sample in samples)
    {
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
