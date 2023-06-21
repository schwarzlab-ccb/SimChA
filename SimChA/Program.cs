using System.Diagnostics;
using System.Globalization;
using CommandLine;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.IO;
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

SimParams simParams;
if (options.ConfigFile != "")
{
    simParams = FileIO.ReadSimParams(options.ConfigFile);
}
else
{
    int seed = new Random().Next();
    var fitness = new FitnessParams(1, 1, 1);
    simParams = new SimParams(seed, SexEnum.Both, 1, Distribution.Uniform, fitness);
}

var rnd = new Random(simParams.Seed);
var files = new FileIO(options.OutputPath);
var genRef = FileIO.GetGenRef(options.DataFolder);
files.WriteSimParams(simParams);

var watch = new Stopwatch();
watch.Start();
List<Sample> samples;
var genContents = options.UseVariants
    ? FileIO.ReadFasta(options.VariantsFile).ToList()
    : new List<GenContents>();
if (execMode == ExecMode.Profiles)
{
    Console.WriteLine("Reading profiles:");
    var profiles = FileIO.ReadProfiles(genRef, options.CNProfiles);
    samples = Simulator.SamplesFromProfiles(profiles);
}
else
{
    var sigs = Validators.ValidateSignatures(simParams.Signatures);
    Console.WriteLine("Computing mutations:");
    var simulator = new Simulator(rnd, geneLists, genContents);
    if (execMode == ExecMode.Tree)
    {
        var inClones = FileIO.ReadClones(options.CloneTreeFile, options.UseMCMC);
        var (cnEventPs, mixture) = Converters.PropagateSigs(sigs);
        string sampleName = Path.GetFileNameWithoutExtension(options.CloneTreeFile);
        var treeSample = new Sample(sampleName, Sampling.GetBinarySex(rnd, simParams.Sex), inClones, cnEventPs, mixture);
        samples = new List<Sample> {treeSample};
    }
    else
    {
        samples = Converters.MakeSamples(rnd, options.Repeats, simParams.EventCount, simParams.Distribution, sigs,
            simParams.Sex);
    }

    foreach (var sample in samples)
    {
        // Monte Carlo sampling of copy-number altering events
        if (options.UseMCMC)
        {
            if (simParams.MCParams == null)
            {
                throw new Exception("Error: MCParams not set. Cannot perform MC sampling. Please set MCParams.");
            }

            simulator = new MCSimulator(rnd, genRef, genContents, simParams.Fitness, simParams.MCParams);
            simulator.SampleEvents(sample);
        }
        else
        {
            simulator.SampleEvents(sample);
        }
    }
}

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

Console.WriteLine("");

try
{
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
}
catch (Exception e)
{
    Console.WriteLine($"Failed with exception {e.Message}. Stack: {e.StackTrace}");
    return e.HResult;
}

watch.Stop();
Console.WriteLine($"Total time: {TimeSpan.FromMilliseconds(watch.ElapsedMilliseconds)}");

return 0;