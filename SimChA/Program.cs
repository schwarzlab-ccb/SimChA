using System.Diagnostics;
using CommandLine;
using SimChA;
using SimChA.Computation;
using SimChA.IO;
using SimChA.Simulation;

var options = Parser.Default.ParseArguments<CmdOptions>(args);
options.WithNotParsed(o =>
{
    Console.WriteLine("Exiting");
    o.ToList().ForEach(Console.Write); // Write out errors
    Environment.Exit(1);
});

SimParams simParams;
if (options.Value.ConfigFile != "")
{
    simParams = FileIo.SimParamsFromFile(options.Value.ConfigFile);
}
else
{
    var aberrations = AberrationsInfo.DefaultAberrations();
    simParams = SimParams.CreateSimParams(new Random().Next(), true, aberrations);
}

string[] newickString = Array.Empty<string>();
if (options.Value.NewickFile != "")
{
    newickString = FileIo.GetStringFromNewick(options.Value.NewickFile);
}

var random = new Random(simParams.Seed);
FileIo files;
try
{
    files = new FileIo(options.Value.OutputPath);
    files.WriteSimParams(simParams);
}
catch (Exception e)
{
    Console.WriteLine($"Failed to write to disk with error: {e.Message}");
    return 2;
}

try
{
    var watch = new Stopwatch();
    watch.Start();

    LcaTreeBuilder.IsNewick = true;
    var clones = Newick.ParseNewickTree(newickString, simParams.IsFemale);
    var aberrations = new AberrationsInfo(simParams);
    Simulator.GetMutationsNewick(clones[0], clones, aberrations, random);
    Console.WriteLine("Mutations generated");
    
    var selectClones = clones.Where(c => c.IsAlive).Shuffle(random).ToList();
    var lcaTree = LcaTreeBuilder.BuildTree(clones, selectClones);

    files.WriteClones(clones);
    files.WriteCopyNumbers(clones);
    files.WriteParentTree(lcaTree);
    
    watch.Stop();
    Console.WriteLine($"Total time: {TimeSpan.FromMilliseconds(watch.ElapsedMilliseconds)}");
}
catch (Exception e)
{
    Console.WriteLine($"Failed with exception {e.Message}. Stack: {e.StackTrace}");
    return e.HResult;
}

return 0;