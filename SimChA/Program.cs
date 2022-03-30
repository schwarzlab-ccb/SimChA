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
        Checkpoints = false,
        // Function
        MultiplicativeFitness = false,
        StochasticCellLife = true,
        FitnessType = FitnessSampleType.Constant,
        Seed = new Random().Next(),
        // Experiment
        PopLimit = 1_000_000_000,
        StepLimit = 100_000,
        CutOff = 0.01f,
        Repeats = 1,
        InitPop = 100,
        // Model
        BirthRate = 0.01,
        MutationRate = 0.00008,
        FitnessMean = 0.05,
        Confinement = 0.1,
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

var globalWatch = new System.Diagnostics.Stopwatch();
globalWatch.Start();

for (int repeatId = 0; repeatId < simParams.Repeats; repeatId++)
{
    var watch = new System.Diagnostics.Stopwatch();
    watch.Start();

    int firstCp = (int) Math.Ceiling(Math.Log2(simParams.InitPop));
    int lastCp = (int) Math.Ceiling(Math.Log2(simParams.PopLimit));
    var checkpoints = simParams.Checkpoints
            ? Enumerable.Range(firstCp, lastCp - firstCp + 1).Select(mag => (int) Math.Pow(2, mag)).ToList()
            : new List<int>();

    Console.WriteLine(string.Join("", Enumerable.Repeat("*", 100)));
    Console.WriteLine($"* Simulation {repeatId + 1}/{simParams.Repeats}");
    Console.WriteLine($"* Sim with seed {simParams.Seed}, genome length  {ReferenceGenome.TotalLength(true)}");

    // Simulation
    var simulator = new Simulator(simParams, random);
    int stepNo = 0;
    var popSizes = new List<(long total, long alive)>();
    popSizes.Add((CellSampling.PopulationSize(simulator.Populations), CellSampling.AliveCount(simulator.Populations)));
    var EndCond = ()
        => !(popSizes.Last().total <= simParams.PopLimit && popSizes.Last().alive > 0 && stepNo < simParams.StepLimit);
    do
    {
        int cloneCount = simulator.FlatPops.Count();
        int popCount = simulator.Populations.Count;
        Console.Write(($"step: {++stepNo:D3}, " +
                       $"populations: {popCount}, " +
                       $"subClones: {cloneCount}, " +
                       $"alive SC: {simulator.AliveSC}, " +
                       $"cells: {popSizes.Last().total:N0}, " +
                       $"alive: {popSizes.Last().alive:N0}").PadRight(160) + "\r");
        simulator.Step();
        popSizes.Add((
            CellSampling.PopulationSize(simulator.Populations),
            CellSampling.AliveCount(simulator.Populations)
        ));

        if (EndCond() || (checkpoints.Any() && popSizes.Last().total > checkpoints.First()))
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
            int stepId = 0;
            if (checkpoints.Any())
            {
                stepId = (int) Math.Log2(checkpoints.First()) - firstCp;
                checkpoints.RemoveAt(0);
            }
            
            var result = new ResultSummary(repeatId, stepId, connectedTree,
                aboveCutOff, cloneCount, sample.Count, stepNo, popSizes, popCount);
            files.AddToSummary(result);

            // Result
            if (EndCond())
            {
                var vaf = TreeAnalysis.ComputeVAF(connectedTree);
                files.WriteFinalOutput(repeatId, sample, lcaTree, aboveCutOff, connectedTree, vaf,
                    popSizes.Last().total, simParams.AliveOnly);
                Console.WriteLine("Final Result:".PadRight(160));
                Console.WriteLine(result.ToText());
            }
        }
    } while (!EndCond());
    
    watch.Stop();
    Console.WriteLine($"Execution Time: {TimeSpan.FromMilliseconds(globalWatch.ElapsedMilliseconds)}\n");
}

files.CopySummary();
globalWatch.Stop();
Console.WriteLine($"Total time: {TimeSpan.FromMilliseconds(globalWatch.ElapsedMilliseconds)}");

return 0;