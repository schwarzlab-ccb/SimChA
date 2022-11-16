using System.Diagnostics;
using CommandLine;
using SimChA;
using SimChA.Computation;
using SimChA.IO;
using SimChA.Simulation;
using SimChA.DataTypes;

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
    simParams = SimParams.CreateSimParams(new Random().Next(), true, defaultAberrs);
}

string newickString = "";
if (options.Value.NewickFile != "")
{
    newickString = FileIO.GetStringFromNewick(options.Value.NewickFile);
}
List<Gen> genes = FileIO.ReadGenes(options.Value.GenesFolder);

Console.WriteLine("Computing mutations.");
var watch = new Stopwatch();
watch.Start();

var random = new Random(simParams.Seed);
var files = new FileIO(options.Value.OutputPath);

var clones = (options.Value.NewickFile != "")
    ? Newick.ParseNewick(newickString, simParams.IsFemale)
    : Simulator.GetClonePair(options.Value.Distance , true);
var aberrationsIfo = new AberrationsInfo(simParams);
Simulator.AssignMutationsRecursive(clones[0], clones, aberrationsIfo, random);
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
    files.WriteTSV();
}
catch (Exception e)
{
    Console.WriteLine($"Failed with exception {e.Message}. Stack: {e.StackTrace}");
    return e.HResult;
}

return 0;