using System.Diagnostics;
using System.Globalization;
using CommandLine;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Simulation;

Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
var watch = new Stopwatch();
watch.Start();

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
    simParams = new SimParams(seed, true, 1, Distribution.Uniform, GenomeAssembly.hg38, fitness, null, null);
}

HGRef.Assembly = simParams.Assembly;
var rnd = new Random(simParams.Seed);
var files = new FileIO(options.Value.OutputPath);
var geneLists = FileIO.ReadGeneLists(options.Value.GenesFolder, simParams.SexXX, HGRef.Assembly);
files.WriteSimParams(simParams);

// Obtain clones
List<Sample> samples;
if (options.Value.CNProfiles != "")
{
    var profiles = FileIO.ReadProfiles(options.Value.CNProfiles);
    samples = Simulator.SamplesFromProfiles(profiles);
}
else
{
    var sigs = Validators.ValidateSignatures(simParams.Signatures);
    Console.WriteLine("Computing mutations.");
    var simulator = new Simulator(rnd, simParams.Fitness, simParams.MCParams, geneLists);
    if (options.Value.NewickFile != "")
    {
        var inClones = FileIO.ReadClones(options.Value.NewickFile, options.Value.UseMCMC);
        var eventPs = Converters.PropagateSigs(sigs);
        var treeSample = new Sample(options.Value.NewickFile, simParams.SexXX, inClones, eventPs);
        samples = new List<Sample> { treeSample };
    }
    else
    {
        samples = Converters.MakeSamples(rnd, options.Value.Repeats, simParams.EventCount, simParams.Distribution, sigs);
    }
    foreach (var sample in samples)
    {
        // Monte Carlo sampling of copy-number altering events
        if (options.Value.UseMCMC && simParams.MCParams == null)
        {
            throw new Exception("Error: MCParams not set. Cannot perform MC sampling. Please set MCParams.");
        }
        simulator.SampleEvents(sample, options.Value.UseMCMC);
    }
}

// Fitness data
Console.WriteLine("Computing fitness.");
var sampleStats = new Dictionary<string, List<CloneStat>>();
int counter = 0;
int cloneCount = samples.Sum(s => s.Clones.Count);
foreach (var sample in samples)
{
    var cloneStatsList = new List<CloneStat>();
    foreach (var clone in sample.Clones)
    {
        Console.Write($"\r{++counter}/{cloneCount} clones processed.");
        var cloneStats = CNProfile.GetCloneStats(clone, geneLists, simParams.Fitness, sample.Kars);
        cloneStatsList.Add(cloneStats);
    }
    sampleStats[sample.SampleId] = cloneStatsList;
}

Console.WriteLine();
try
{
    files.WriteSamples(samples);
    files.WriteClones(samples);
    files.WriteCopyNumbers(samples);
    files.WriteFitness(sampleStats);
    // Check if cnEvents was assigned
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