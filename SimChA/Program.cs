using System.Diagnostics;
using CommandLine;
using SimChA.Computation;
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
    const bool isFemale = true;
    const float pStress = 000_01f;
    const float pTsgOg = 000_1f;
    const float pEssential = 000_01f;
    var pAberrs = AberrationsInfo.DefaultAberrations();
    simParams = SimParams.CreateSimParams(new Random().Next(), isFemale, pStress, pTsgOg, pEssential, pAberrs);
}

var files = new FileIO(options.Value.OutputPath);
var rnd = new Random(simParams.Seed);

var geneLists = FileIO.ReadGeneLists(options.Value.GenesFolder, simParams.IsFemale);

string newickString = "";
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
    var result = new Dictionary<string, double>();
    int counter = 0;
    foreach ((string sample, var kar) in cnas)
    {
        Console.Write("\r" + counter++ + " / " + cnas.Count + " samples processed.");
        double fitness = Fitness.Calculate(kar, geneLists, simParams);
        result.Add(sample, fitness);
    }
    Console.WriteLine("Writing to disk.");
    files.WriteSampleFitness(result);
}
else
{
    Console.WriteLine("Computing mutations.");

    var clones = options.Value.NewickFile != ""
        ? Parsers.ParseNewick(newickString, simParams.IsFemale)
        : Simulator.MakeClonePair(options.Value.Distance, true);
    var aberrationsInfo = new AberrationsInfo(simParams);
    var simulator = new Simulator(aberrationsInfo, rnd, simParams, geneLists);
    var aberrations = simulator.AssignMutations(clones[0], clones);

    // TODO: do not remove the diploid clone if a newick file is provided
    var selectClones = clones.Where(c => c.CloneId != 0).ToList();

    Console.WriteLine("Writing to disk.");
    try
    {
        files.WriteClones(selectClones);
        files.WriteCopyNumbers(selectClones, simParams.IsFemale);
        files.WriteSimParams(simParams);
        files.WriteTSV(aberrations);
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