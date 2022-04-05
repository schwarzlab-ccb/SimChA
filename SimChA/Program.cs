using SimChA.DataTypes;
using SimChA.Simulation;
using CommandLine;
using SimChA.Computation;
using SimChA.IO;

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
        AliveOnly = true,
        Checkpoints = true,
        // Function
        MultiplicativeFitness = false,
        StochasticCellLife = true,
        FitnessType = FitnessSampleType.Constant,
        Seed = new Random().Next(),
        // Experiment
        PopLimit = 100_000_000,
        StepLimit = 100_000,
        CutOff = 0.01f,
        Repeats = 1,
        InitPop = 1000,

        // Model
        BirthRate = 0.01,
        MutationRate = 0.00016,
        FitnessMean = 0.5,
        Confinement = 0.0,
        InitMut = 0,
    };
}

var random = new Random(simParams.Seed);
FileIO files;
bool isRepeated = simParams.Repeats > 1;
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
    var globalWatch = new System.Diagnostics.Stopwatch();
    globalWatch.Start();
    
    int tryOut = 0;
    for (int repeatId = 0; repeatId < simParams.Repeats; repeatId++)
    {
        var watch = new System.Diagnostics.Stopwatch();
        watch.Start();

        // Simulation
        int stepNo = 0;
        int checkpointId = 0;
        var simulator = new Simulator(simParams, random);
        var checkpoints = Utility.CreateCheckpoints(simParams);
        var popSizes = new List<(long total, long alive)> { CellSampling.PopState(simulator.FlatPops.ToList()) };
        var EndCondFunc = () =>
            !(popSizes.Last().total <= simParams.PopLimit
              && popSizes.Last().alive > 0
              && stepNo < simParams.StepLimit);
        do
        {
            simulator.Step();
            
            var flatPops = simulator.FlatPops.ToList();
            int popCount = simulator.Populations.Count;
            popSizes.Add(CellSampling.PopState(flatPops));
            if (popSizes.Last().alive <= 0 && popSizes.Last().total < simParams.InitPop)
            {
                break;
            }
            Console.Write(($"Sim: {repeatId + 1}.{tryOut}/{simParams.Repeats}, "+
                           $"step: {++stepNo:D3}, " +
                           $"populations: {popCount}, " +
                           $"subClones: {flatPops.Count}, " +
                           $"alive SC: {simulator.AliveSC}, " +
                           $"cells: {popSizes.Last().total:N0}, " +
                           $"alive: {popSizes.Last().alive:N0}").PadRight(160) +
                          (options.Value.Newline ? "\n" : "\r"));

            if (EndCondFunc() || checkpoints.Any() && popSizes.Last().total > checkpoints.First())
            {
                // Analysis
                var cutOff = popSizes.Select(pair => (long)Math.Ceiling(pair.alive * simParams.CutOff)).ToList();
                var aboveCutOff = simulator.FlatPops.Where(sc
                    => Enumerable.Range(0, popSizes.Count).Any(g => cutOff[g] <= sc.AliveAtGen(g))).ToList();
                var lcaTree = LCATreeBuilder.Builtree(simulator.FlatPops, aboveCutOff);
                var connectedTree = ConnectedTreeBuilder.BuildTree(simulator.FlatPops, aboveCutOff);
                var treeNodes = lcaTree.Nodes.Select(n => n.Id).ToList();
                var sample = simulator.FlatPops.Where(sc => treeNodes.Contains(sc.CloneId)).ToList();

                // Summary
                if (checkpoints.Any())
                {
                    checkpointId++;
                    checkpoints.RemoveAt(0);
                }

                var result = new ResultSummary(repeatId, checkpointId, connectedTree,
                    aboveCutOff, flatPops.Count, simulator.AliveSC, sample.Count, stepNo, popSizes, popCount);
                files.AddToSummary(result);

                // Result
                if (EndCondFunc())
                {
                    var vaf = TreeAnalysis.ComputeVAF(connectedTree);
                    files.WriteFinalOutput(repeatId, sample, lcaTree, aboveCutOff, connectedTree, vaf,
                        popSizes.Last().total, simParams.AliveOnly);
                    Console.WriteLine($"Sim: {repeatId + 1}.{tryOut}/{simParams.Repeats} result:".PadRight(160));
                    Console.WriteLine(result.ToText());
                }
            }
        } while (!EndCondFunc());

        if (popSizes.Last().alive <= 0 && popSizes.Last().total < simParams.InitPop)
        {
            repeatId--;
            tryOut++;
            continue;
        }

        watch.Stop();
        Console.WriteLine($"Execution Time: {TimeSpan.FromMilliseconds(globalWatch.ElapsedMilliseconds)}");
        Console.WriteLine(string.Join("", Enumerable.Repeat("*", 100)));
        tryOut = 0;
    }

    files.CopySummary();
    globalWatch.Stop();
    
    Console.WriteLine($"Total time: {TimeSpan.FromMilliseconds(globalWatch.ElapsedMilliseconds)}");
}
catch (Exception e)
{
    Console.WriteLine($"Failed with exception {e.Message}");
}

return 0;