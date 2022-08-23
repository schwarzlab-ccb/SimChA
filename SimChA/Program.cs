using System.Diagnostics;
using CommandLine;
using SimChA;
using SimChA.Computation;
using SimChA.DataTypes;
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
    simParams = FileIO.SimParamsFromFile(options.Value.ConfigFile);
}
else
{
    simParams = new SimParams
    {
        Seed = new Random().Next(),
        // Experiment
        CloneTarget = 10_000,
        SelectionSize = 25,
        MaxSteps = 100_000,
        Reps = 1,

        // Model
        IsFemale = true,
        Turnover = 0.01,
        MutationProb = 0.25,
        DeathRate = 0.95,

        // Initialization
        StartMut = 1,
        StartPop = 1,
        AberrationRates =
        {
            [AberrationEnum.InternalDeletion] = 50f,
            [AberrationEnum.InternalDuplication] = 50f,
            [AberrationEnum.Translocation] = 20f,
            [AberrationEnum.TailDeletion] = 15f,
            [AberrationEnum.BreakageFusionBridge] = 10f,
            [AberrationEnum.Inversion] = 10f,
            [AberrationEnum.Missegregation] = 5f,
            [AberrationEnum.Duplication] = 5f,
            [AberrationEnum.Chromothripsis] = 1f,
            [AberrationEnum.WholeGenomeDoubling] = 1f
        }
    };
}

string[] newickString = {};
if(options.Value.NewickFile != ""){
    newickString = FileIO.GetStringFromNewick(options.Value.NewickFile);
}

var random = new Random(simParams.Seed);
FileIO files;
bool isRepeated = simParams.Reps > 1;
try
{
    files = new FileIO(options.Value.OutputPath, isRepeated);
    files.WriteSimParams(simParams);
}
catch (Exception e)
{
    Console.WriteLine($"Failed to write to disk with error: {e.Message}");
    return 2;
}

try
{
    var globalWatch = new Stopwatch();
    globalWatch.Start();
    if(options.Value.NewickFile == ""){
        int tryNo = 0;
        for (int repeatId = 0; repeatId < simParams.Reps; repeatId++)
        {
            var watch = new Stopwatch();
            watch.Start();
            // Simulation

            string lastLine = "";
            var simulator = new Simulator(simParams, random);
            do
            {
                simulator.StepTree();
                int lastSize = lastLine.Length; // Only pad what you need
                double prog = (double)simulator.Clones.Count / simParams.CloneTarget;
                lastLine = $"sim: {repeatId}.{tryNo}/{simParams.Reps}, " +
                        $"step: {simulator.StepNo:D3}, " +
                        $"prog: {prog:P}, " +
                        $"SC_total: {simulator.Clones.Count}, " +
                        $"SC_alive: {simulator.AliveClones},";
                Console.Write(lastLine.PadRight(lastSize) + (options.Value.Newline ? "\n" : "\r"));
            } while (simulator.Clones.Count < simParams.CloneTarget && simulator.StepNo < simParams.MaxSteps && simulator.AliveClones > 0);

            if (simulator.AliveClones <= 0)
            {
                repeatId--;
                tryNo++;
                continue;
            }

            Console.WriteLine($"Starting Mutations");
            //simulator.GetMutations(simulator.Clones[0]); //Starting with root

            // Analysis
            var selectClones = simulator.Clones
                .Where(c => c.IsAlive)
                .Shuffle(random)
                .Take(simParams.SelectionSize)
                .ToList();
            var lcaTree = LCATreeBuilder.Builtree(simulator.Clones, selectClones);
            var treeNodes = lcaTree.Nodes.Select(n => n.Id).ToList();
            var sample = simulator.Clones.Where(sc => treeNodes.Contains(sc.CloneId)).ToList();
            var snps = SNPBuilder.CreateSNPs(random, simParams.IsFemale, 100); // snps are shared between all subclones and therefore are created only once

            Console.Write("".PadRight(lastLine.Length) + "\r");
            files.WriteClones(sample);
            files.WriteParentTree(lcaTree);
            files.StoreCopy(repeatId);
            files.WriteClones(sample);
            files.WriteRawData(random, sample, snps, simParams.IsFemale);
            files.WriteCopyNumbers(sample);
        }

        files.RetrieveConfig();
        globalWatch.Stop();

        Console.WriteLine($"Total time: {TimeSpan.FromMilliseconds(globalWatch.ElapsedMilliseconds)}");
    }
    else{
        var watch = new Stopwatch();
        LCATreeBuilder.isNewick = true;
        watch.Start();
        var simulator = new Simulator(simParams, random);
        simulator.BuildCloneFromNewick(newickString);
        simulator.GetMutationsNewick(simulator.Clones[0]);
        Console.WriteLine("Mutations generated");
        var selectClones = simulator.Clones
                .Where(c => c.IsAlive)
                .Shuffle(random)
                .Take(simParams.SelectionSize)
                .ToList();        
        var lcaTree = LCATreeBuilder.Builtree(simulator.Clones, selectClones);
        var treeNodes = lcaTree.Nodes.Select(n => n.Id).ToList();
        var sample = simulator.Clones.Where(sc => treeNodes.Contains(sc.CloneId)).ToList();


        files.WriteClones(simulator.Clones);
        files.WriteCopyNumbers(simulator.Clones);
        files.WriteParentTree(lcaTree);
        watch.Stop();
        Console.WriteLine($"Total time: {TimeSpan.FromMilliseconds(globalWatch.ElapsedMilliseconds)}");
    }
}
catch (Exception e)
{
    Console.WriteLine($"Failed with exception {e.Message}. Stack: {e.StackTrace}");
    return e.HResult;
}

return 0;