using System.Diagnostics;
using CommandLine;
using SimChA;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Simulation;

ParserResult<CmdOptions>? options = Parser.Default.ParseArguments<CmdOptions>(args);
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
    simParams = new SimParams
    {
        Seed = new Random().Next(),
        IsFemale = true,

        AberrationRates =
        {
            [AberrationEnum.InternalDeletion] = 50f,
            [AberrationEnum.InternalDuplication] = 50f,
            [AberrationEnum.Translocation] = 20f,
            [AberrationEnum.TailDeletion] = 15f,
            [AberrationEnum.BreakageFusionBridge] = 10f,
            [AberrationEnum.Inversion] = 10f,
            [AberrationEnum.ChromDeletion] = 5f,
            [AberrationEnum.ChromDuplication] = 5f,
            [AberrationEnum.Chromothripsis] = 1f,
            [AberrationEnum.WholeGenomeDoubling] = 1f
        }
    };
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
    var simulator = new Simulator(simParams, random);
    simulator.BuildCloneFromNewick(newickString);
    simulator.GetMutationsNewick(simulator.Clones[0]);
    Console.WriteLine("Mutations generated");
    var selectClones = simulator.Clones.Where(c => c.IsAlive).Shuffle(random).ToList();
    ParentTree lcaTree = LcaTreeBuilder.BuildTree(simulator.Clones, selectClones);

    files.WriteClones(simulator.Clones);
    files.WriteCopyNumbers(simulator.Clones);
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