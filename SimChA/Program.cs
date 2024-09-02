using System.Diagnostics;
using System.Globalization;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Optimization;
using SimChA.Simulation;
using CommandLine;

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

var simParams = FileIO.ReadSimParams(options.ConfigFile);

var rnd = new Random(simParams.Seed);
var files = new FileIO(options.OutputPath);
bool parseGenContents = execMode == ExecMode.ParseGenContents;
bool includeSexChromosomes = !options.AutosomesOnly;
var genRef = FileIO.GetGenRef(options.DataFolder, includeSexChromosomes, parseGenContents);

var watch = new Stopwatch();
watch.Start();
List<Sample> samples;
var simulator = new Simulator(rnd, genRef);

switch (execMode)
{
    case ExecMode.BinSamples:
    {
        string cnProfiles = options.CNProfiles ?? throw new Exception("Error: No CN profiles provided. Cannot bin samples.");
        Console.WriteLine("Reading observed data:");
        var profiles = FileIO.ReadProfiles(genRef, cnProfiles);
        var observedSamples = Simulator.SamplesFromProfiles(profiles);
        var binner = new Binner(genRef);
        Console.WriteLine("Binning samples:");
        var binnedSamples = binner.GetBinnedCNProfiles(observedSamples);
        Console.WriteLine("Writing binned samples:");
        files.WriteCopyNumbers(binnedSamples);
        return 0;
    }
    case ExecMode.RunOptimization:
    {
        var optParams = simParams.OptimizationParams ?? throw new Exception("Error: OptimizationParams not set. Cannot perform optimization. Please set OptimizationParams.");
        if (optParams.Mode != "Events" && optParams.Mode != "Fitness")
        {
            throw new Exception("Error: OptimizationParams.Mode not recognized. Please set OptimizationParams.Mode to either Events or Fitness.");
        }
    
        var optimizer = simParams.OptimizationParams.Mode == "Events" 
            ? new Optimizer(simParams, rnd, options.Repeats, genRef, includeSexChromosomes, files)
            : new FitnessOptimizer(simParams, rnd, options.Repeats, genRef, includeSexChromosomes, files);

        if (options.TargetParams != "")
        {
            Console.WriteLine("Generating synthetic data");
            var targetParams = FileIO.ReadSimParams(options.TargetParams);
            optimizer.InitializeObservations(targetParams);
        }
        else if (options.CNProfiles != "" && options.EventCounts != "")
        {
            Console.WriteLine("Reading observed data:");
            var profiles = FileIO.ReadProfiles(genRef, options.CNProfiles);
            var observedSamples = Simulator.SamplesFromProfiles(profiles);
            var eventCounts = FileIO.ReadEventCounts(options.EventCounts);
            optimizer.InitializeObservations(observedSamples, eventCounts);
        }
        else 
        {
            throw new Exception("Error: No target parameters (synthetic data) or bootstrap file (observed data) provided. Cannot perform event optimization without a target data set.");
        }
    
        // Perform the optimization
        var outParams = optimizer.Optimize();
        files.WriteSimParams(outParams);
        watch.Stop();
        Console.WriteLine($"Total time: {TimeSpan.FromMilliseconds(watch.ElapsedMilliseconds)}");
        Console.WriteLine("Optimization finished");
        return 0;
    }
    case ExecMode.Profiles:
    {
        Console.WriteLine("Reading profiles:");
        var profiles = FileIO.ReadProfiles(genRef, options.CNProfiles);
        samples = Simulator.SamplesFromProfiles(profiles);
        if (options.FitnessLandscape)
        {
            Console.WriteLine("Computing fitness landscape:");
            FitnessLandscape.GenerateFitnessLandscape(genRef, simParams, samples, files);
        }
        break;
    }
    
    case ExecMode.None:
    case ExecMode.Tree:
    case ExecMode.Repeats:
    case ExecMode.ParseGenContents:
    case ExecMode.Bootstrap:
    case ExecMode.UseMCMC:
    default:
    {
        var sigs = simParams.Signatures ?? throw new Exception("Error: Signatures not set. Cannot perform simulation without signatures. Please set Signatures in the config file.");
        Validators.ValidateSignatures(simParams.Signatures);
        Console.WriteLine("Computing mutations:");
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
            var mcParams = simParams.MCParams ?? throw new Exception("Error: MCParams not set. Cannot perform MCMC sampling. Please set MCParams in the config file.");
            if (options.CNProfiles != "" && options.EventCounts != "")
            {
                var profiles = FileIO.ReadProfiles(genRef, options.CNProfiles);
                var eventCounts = FileIO.ReadEventCounts(options.EventCounts);
                var fitnessList = simulator.FitnessListFromSamples(simParams, profiles, eventCounts);
                samples = Converters.MakeSamples(rnd, options.Repeats, simParams.EventCount, simParams.EventDist, sigs, simParams.Sex, fitnessList);
            }
            else
            {
                var mcTarget =  simParams.MCTarget ?? throw new Exception("Error: MCTarget not set. Cannot perform MC sampling. Please set MCTarget in the config file.");
                samples = Converters.MakeSamples(rnd, options.Repeats, simParams.EventCount, simParams.EventDist, sigs, simParams.Sex, mcTarget);
            }
            simulator = new MCSimulator(rnd, genRef, simParams.Fitness, mcParams);
        }
        else
        {
            var mcTarget =  simParams.MCTarget ?? throw new Exception("Error: MCTarget not set. Cannot perform MC sampling. Please set MCTarget in the config file.");
            var eventCounts = options.EventCounts != "" ? FileIO.ReadEventCounts(options.EventCounts) : new Dictionary<string, int>();
            samples = Converters.MakeSamples(rnd, options.Repeats, simParams.EventCount, simParams.EventDist, sigs, simParams.Sex, mcTarget, eventCounts);
        }

        foreach (var sample in samples)
        {
            simulator.SampleEvents(sample);
        }
        break;
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

// Write output
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
