using System.Diagnostics;
using CommandLine;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Simulation;

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

var cnas = new Dictionary<string, Karyotype>();
if (options.Value.CNProfiles != "")
{
    cnas = FileIO.ReadCopyNumbers(options.Value.CNProfiles);
}

var watch = new Stopwatch();
watch.Start();

if (options.Value.CNProfiles != "")
{
    Console.WriteLine("Computing fitness.");
    var results = new List<ProfileStats>();
    int counter = 0;
    foreach ((string sample, var kar) in cnas)
    {
        Console.Write($"\r{counter++}/{cnas.Count} samples processed.");
        var profileStats = CNProfile.GetProfileStats(sample, kar, geneLists, simParams.Fitness);
        results.Add(profileStats);
    }
    Console.WriteLine("Writing to disk.".PadRight(80));
    files.WriteSampleFitness(results);
}
else
{
    Parsers.ValidateSignatures(simParams.Signatures);
    Console.WriteLine("Computing mutations.");

    var clones = options.Value.NewickFile != ""
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
watch.Stop();
Console.WriteLine($"Total time: {TimeSpan.FromMilliseconds(watch.ElapsedMilliseconds)}");

return 0;