using System.Diagnostics;
using System.Globalization;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.IO;
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
var fitParams = simParams.Fitness ?? throw new Exception("Error: Fitness parameters are missing. Please set Fitness in the config file.");

var rnd = new Random(simParams.Seed);
var files = new FileIO(options.OutputPath);
var genRef = FileIO.GetGenRef(options.DataFolder, options.ShouldParseGenome);

var watch = new Stopwatch();
watch.Start();
List<Sample> samples;
Simulator simulator;
if (options.UseMCMC)
{
    var mcParams = simParams.MCParams ?? throw new Exception("Error: MCParams not set. Cannot perform MC sampling. Please set MCParams in the config file.");
    simulator = new MCSimulator(rnd, genRef, fitParams, mcParams, files);
}
else
{
    simulator = new Simulator(rnd, genRef);
}

switch (execMode)
{
    case ExecMode.Profiles:
    {
        Console.WriteLine("Reading profiles:");
        var profiles = FileIO.ReadProfiles(genRef, options.CNProfiles, simParams.AutosomesOnly);
        // TODO: Needs to implement autosomes only
        samples = Simulator.SamplesFromProfiles(profiles);
        break;
    }
    
    case ExecMode.Tree:
        var treeSigs = simParams.Signatures ?? throw new Exception("Error: Signatures not set. Cannot perform simulation without signatures. Please set Signatures in the config file.");
        Validators.ValidateSignatures(treeSigs);
        Console.WriteLine("Computing mutations for tree:");
        var inClones = FileIO.ReadClones(options.CloneTreeFile, options.UseMCMC);
        var (cnEventPs, mixture) = Converters.PropagateSigs(treeSigs);
        string sampleName = Path.GetFileNameWithoutExtension(options.CloneTreeFile);
        var sex = simParams.AutosomesOnly ? SexEnum.None : Sampling.GetSex(rnd, simParams.Sex);
        var treeSample = new Sample(sampleName, sex, inClones, cnEventPs, mixture);
        samples = new List<Sample> {treeSample};
        samples.ForEach(simulator.SampleEvents);
        break;
    
    case ExecMode.Repeats:
    default:
        var repSigs = simParams.Signatures ?? throw new Exception("Error: Signatures not set. Cannot perform simulation without signatures. Please set Signatures in the config file.");
        Validators.ValidateSignatures(repSigs);
        Console.WriteLine("Computing mutations for individual samples:");
        if (options.UseMCMC)
        {
            if (options.CNProfiles != "" && options.EventCounts != "")
            {
                var profiles = FileIO.ReadProfiles(genRef, options.CNProfiles, simParams.AutosomesOnly);
                var eventCounts = FileIO.ReadEventCounts(options.EventCounts);
                var fitnessList = simulator.FitnessListFromSamples(simParams, profiles, eventCounts);
                samples = Converters.MakeSamples(rnd, options.Repeats, simParams.EventCount, simParams.EventDist, repSigs, simParams.Sex, simParams.AutosomesOnly, fitnessList);
            }
            else if (options.EventCounts != "" && simParams.MCParams != null && !simParams.MCParams.MatchFitness)
            {
                var eventCounts = FileIO.ReadEventCounts(options.EventCounts);
                samples = Converters.MakeSamples(rnd, options.Repeats, simParams.EventCount, simParams.EventDist, repSigs, simParams.Sex, simParams.AutosomesOnly);
            }
            else
            {
                var mcTarget =  simParams.MCTarget ?? throw new Exception("Error: MCTarget not set. Cannot perform MC sampling. Please set MCTarget in the config file.");
                samples = Converters.MakeSamples(rnd, options.Repeats, simParams.EventCount, simParams.EventDist, repSigs, simParams.Sex, simParams.AutosomesOnly, mcTarget);
            }
        }
        else
        {
            var mcTarget =  simParams.MCTarget ?? throw new Exception("Error: MCTarget not set. Cannot perform MC sampling. Please set MCTarget in the config file.");
            var eventCounts = options.EventCounts != "" ? FileIO.ReadEventCounts(options.EventCounts) : new Dictionary<string, int>();
            samples = Converters.MakeSamples(rnd, options.Repeats, simParams.EventCount, simParams.EventDist, repSigs, simParams.Sex, simParams.AutosomesOnly, mcTarget, eventCounts);
        }

        samples.ForEach(simulator.SampleEvents);
        break;
    
    case ExecMode.None:
        throw new Exception("Error: No execution mode set.");
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
        sample.Stats[clone.CloneId] = CNProfile.GetCloneStats(sample, clone, genRef, fitParams, sample.Kars);
    }
}

// TODO split generation of fitness landscape and write to file
if (options.FitnessLandscape)
{
    Console.WriteLine("Computing fitness landscape:");
    FitnessLandscape.GenerateFitnessLandscape(genRef, simParams, samples, files);
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
