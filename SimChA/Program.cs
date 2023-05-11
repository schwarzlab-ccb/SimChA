using System.Diagnostics;
using System.Globalization;
using CommandLine;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Simulation;

Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
var options = Parser.Default.ParseArguments<CmdOptions>(args);
options.WithNotParsed(o =>
{
    Console.WriteLine("Exiting");
    o.ToList().ForEach(Console.WriteLine); // Write out errors
    Environment.Exit(1);
});

SimParams simParams;
string configFile = options.Value.ConfigFile;
if (configFile != "")
{
    simParams = FileIO.ReadSimParams(configFile);
}
else
{
    int seed = new Random().Next();
    var fitness = new FitnessParams(1, 1, 1);
    simParams = new SimParams(seed, true, 1, Distribution.Uniform, GenomeAssembly.hg38, fitness);
}
HGRef.Assembly = simParams.Assembly;
var rnd = new Random(simParams.Seed);
var files = new FileIO(options.Value.OutputPath);
var geneLists = FileIO.ReadGeneLists(options.Value.GenesFolder, simParams.SexXX, HGRef.Assembly);
files.WriteSimParams(simParams);

var watch = new Stopwatch();
watch.Start();
List<Sample> samples;
if (options.Value.CNProfiles != "")
{
    Console.WriteLine("Reading profiles:");
    var profiles = FileIO.ReadProfiles(options.Value.CNProfiles);
    samples = Simulator.SamplesFromProfiles(profiles);
}
else
{
    var sigs = Validators.ValidateSignatures(simParams.Signatures);
    Console.WriteLine("Computing mutations:");
    var simulator = new Simulator(rnd, simParams.Fitness, geneLists);
    if (options.Value.CloneTreeFile != "")
    {
        var inClones = FileIO.ReadClones(options.Value.CloneTreeFile, options.Value.UseMCMC);
        var eventPs = Converters.PropagateSigs(sigs);
        string sampleName = Path.GetFileNameWithoutExtension(options.Value.CloneTreeFile);
        var treeSample = new Sample(sampleName, simParams.SexXX, inClones, eventPs);
        samples = new List<Sample> { treeSample };
    }
    else
    {
        samples = Converters.MakeSamples(rnd, options.Value.Repeats, simParams.EventCount, simParams.Distribution, sigs);
    }
    foreach (var sample in samples)
    {
        // Monte Carlo sampling of copy-number altering events
        if (options.Value.UseMCMC)
        {
            if (simParams.MCParams == null)
            {
                throw new Exception("Error: MCParams not set. Cannot perform MC sampling. Please set MCParams.");
            }
            simulator.SampleEvents(sample, simParams.MCParams);
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
int counter = 1;
int total = samples.Sum(s => s.Clones.Count);
foreach (var sample in samples)
{
    foreach (var clone in sample.Clones)
    {
        Console.Write($"\rClone {counter++}/{total}.");
        sample.Stats[clone.CloneId] = CNProfile.GetCloneStats(clone, geneLists, simParams.Fitness, sample.Kars);
    }
}
Console.WriteLine("");

try
{
    files.WriteSamples(samples);
    files.WriteCopyNumbers(samples);
    files.WriteFitness(samples);
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