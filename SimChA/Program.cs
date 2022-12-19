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
    simParams = FileIO.SimParamsFromFile(options.Value.ConfigFile);
}
else
{
    var defaultAberrs = AberrationsInfo.DefaultAberrations();
    simParams = SimParams.CreateSimParams(new Random().Next(), true, 0.00001f, 0.0001f, 0.00001f, defaultAberrs);
}

string newickString = "";
if (options.Value.NewickFile != "")
{
    newickString = FileIO.GetStringFromNewick(options.Value.NewickFile);
}
var (tsgOgList, essentialsList) = FileIO.ReadGenes(options.Value.GenesFolder, simParams.IsFemale);

Console.WriteLine("Computing mutations.");
var watch = new Stopwatch();
watch.Start();

var rnd = new Random(simParams.Seed);
var files = new FileIO(options.Value.OutputPath);

var clones = options.Value.NewickFile != ""
    ? Newick.ParseNewick(newickString, simParams.IsFemale)
    : Simulator.MakeClonePair(options.Value.Distance , true);
var aberrationsInfo = new AberrationsInfo(simParams);
var simulator = new Simulator(aberrationsInfo, rnd, simParams, tsgOgList, essentialsList);
var aberrations = simulator.AssignMutations(clones[0], clones);
var lcaTree = LcaTreeBuilder.BuildTree(clones);
    
watch.Stop();
Console.WriteLine($"Total time: {TimeSpan.FromMilliseconds(watch.ElapsedMilliseconds)}");

Console.WriteLine("Writing to disk.");
try
{
    files.WriteClones(clones);
    files.WriteCopyNumbers(clones);
    files.WriteParentTree(lcaTree);
    files.WriteSimParams(simParams);
    files.WriteNewickFile(clones);
    files.WriteTSV(aberrations);
}
catch (Exception e)
{
    Console.WriteLine($"Failed with exception {e.Message}. Stack: {e.StackTrace}");
    return e.HResult;
}

return 0;