using System.Diagnostics;
using CommandLine;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Simulation;

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
if (options.Value.ConfigFile != "")
{
    simParams = FileIO.ReadSimParams(options.Value.ConfigFile);
}
else
{
    int seed = new Random().Next();
    var fitness = new FitnessParams(1, 1, 1);
    simParams = new SimParams(seed, true, fitness, null);
}

var rnd = new Random(simParams.Seed);
var files = new FileIO(options.Value.OutputPath);
files.WriteSimParams(simParams);

var geneLists = FileIO.ReadGeneLists(options.Value.GenesFolder, simParams.SexXX);

var newickString = "";
if (options.Value.NewickFile != "")
{
    newickString = FileIO.ReadNewick(options.Value.NewickFile);
}

Dictionary<string, Karyotype> profiles = new();
if (options.Value.CNProfiles != "")
{
    profiles = FileIO.ReadProfiles(options.Value.CNProfiles);
}

List<Clone> clones;
if (options.Value.CNProfiles != "")
{
    clones = Simulator.ClonesFromProfiles(profiles);
}
else
{
    Parsers.ValidateSignatures(simParams.Signatures);
    Console.WriteLine("Computing mutations.");
    
    clones = options.Value.NewickFile != ""
        ? Parsers.ParseNewick(newickString, simParams.SexXX)
        : Simulator.MakeClonePair(options.Value.Distance, simParams.SexXX);
    var simulator = new Simulator(rnd, simParams, geneLists);
    var aberrations = simulator.ApplyEvents(clones[0], clones);
    // TODO: do not remove the diploid clone if a newick file is provided
    var selectClones = clones.Where(c => c.CloneId != 0).ToList();

    Console.WriteLine("Writing to disk.".PadRight(80));
    try
    {
        files.WriteClones(selectClones);
        files.WriteCopyNumbers(selectClones, simParams.SexXX);
        files.WriteEvents(aberrations);
    }
    catch (Exception e)
    {
        Console.WriteLine($"Failed with exception {e.Message}. Stack: {e.StackTrace}");
        return e.HResult;
    }
}

// Fitness data
Console.WriteLine("Computing fitness.");
var results = new List<ProfileStats>();
int counter = 0;
foreach (var clone in clones)
{
    Console.Write($"\r{++counter}/{clones.Count} samples processed.");
    var profileStats = CNProfile.GetProfileStats(clone, geneLists, simParams.Fitness);
    results.Add(profileStats);
}
Console.WriteLine("Writing to disk.".PadRight(80));
files.WriteSampleFitness(results);

watch.Stop();
Console.WriteLine($"Total time: {TimeSpan.FromMilliseconds(watch.ElapsedMilliseconds)}");

return 0;