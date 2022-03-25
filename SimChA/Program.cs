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

var programConfig = new ProgramConfig
{
    MultiplicativeFitness = false,
    StochasticCellLife = false
};

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
        PopLimit = 1_000_000_000,
        StepLimit = 100_000,
        CutOff = 0.01f,
        Repeats = 1,
        DivisionRate = 0.01f,
        MutationRate = 0.00002f,
        FitnessLambdaInv = 0.1f,
        Confinement = 0.05f,
        SplitRate = 0.0f,
        InitialPop = 1,
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

for (int i = 0; i < simParams.Repeats; i++)
{
    var watch = new System.Diagnostics.Stopwatch();
    watch.Start();

    Console.WriteLine(string.Join("", Enumerable.Repeat("*", 100)));
    Console.WriteLine($"* Simulation {i + 1}/{simParams.Repeats}");
    Console.WriteLine($"* Sim with seed {simParams.Seed}, genome length  {ReferenceGenome.TotalLength(true)}");

    // Simulation
    var simulator = new Simulator(simParams, programConfig, random);
    int stepNo = 0;
    var popSizes = new List<(long total, long alive)>();
    popSizes.Add((CellSampling.PopulationSize(simulator.Populations), CellSampling.AliveCount(simulator.Populations)));
    int cloneCount;
    do
    {
        cloneCount = simulator.FlatPops.Count();
        Console.Write(($"Step: {++stepNo:D3}, " +
                       $"populations: {simulator.Populations.Count}, " +
                       $"subClones: {cloneCount}, " +
                       $"alive SC: {simulator.AliveSC}, " +
                       $"cells: {popSizes.Last().total:N0}, " +
                       $"alive: {popSizes.Last().alive:N0}").PadRight(160) + "\r");
        simulator.Step();
        popSizes.Add((
                CellSampling.PopulationSize(simulator.Populations),
                CellSampling.AliveCount(simulator.Populations)
            ));
    } while (popSizes.Last().total <= simParams.PopLimit && popSizes.Last().alive > 0 && stepNo < simParams.StepLimit);

    // Analysis
    var cutOff = popSizes.Select(pair => (long)Math.Ceiling(pair.alive * simParams.CutOff)).ToList();
    var aboveCutOff = simulator.FlatPops.Where(sc
        => Enumerable.Range(0, popSizes.Count).Any(g => cutOff[g] <= sc.AliveAtGen(g))).ToList();
    var lcaTree = LCATreeBuilder.Builtree(simulator.FlatPops, aboveCutOff);
    var connectedTree = ConnectedTreeBuilder.BuildTree(simulator.FlatPops, aboveCutOff);
    // var connectedFullTree = ConnectedTreeBuilder.BuildTree(simulator.FlatPops, simulator.FlatPops.ToList());
    var treeNodes = lcaTree.Nodes.Select(n => n.Id).ToList();
    var sample = simulator.FlatPops.Where(sc => treeNodes.Contains(sc.CloneId)).ToList();
    var vaf = TreeAnalysis.ComputeVAF(connectedTree);
    // var snps = SNPBuilder.CreateSNPs(random, simParams.IsFemale, 100); // snps are shared between all subclones and therefore are created only once

    // Summary
    ResultSummary resultSummary = new();
    (resultSummary.NodeCount, resultSummary.LeafCount, resultSummary.TreeDepth, resultSummary.Branching)
        = TreeAnalysis.ComputeTreeSize(connectedTree);
    // resultSummary.treeBalance = TreeAnalysis.ComputeTreeBalance(connectedFullTree);
    resultSummary.treeBalanceFiltered = TreeAnalysis.ComputeTreeBalance(connectedTree);
    // resultSummary.clonalDiversity = TreeAnalysis.ComputeClonalDiversity(simulator.FlatPops.ToList());
    resultSummary.clonalDiversityFiltered = TreeAnalysis.ComputeClonalDiversity(aboveCutOff);
    // resultSummary.meanDriversPerCell = TreeAnalysis.ComputeMeanDriversPerCell(simulator.FlatPops.ToList());
    resultSummary.meanDriversPerCellFiltered = TreeAnalysis.ComputeMeanDriversPerCell(aboveCutOff);
    resultSummary.SubcloneTotal = cloneCount;
    resultSummary.SubcloneSelect = sample.Count;
    resultSummary.Generations = stepNo;
    resultSummary.AliveCount = popSizes.Last().alive;
    resultSummary.TotalCount = popSizes.Last().total;
    Console.WriteLine("Result:".PadRight(160));
    Console.WriteLine(resultSummary.ToLine());

    // Output
    Console.WriteLine("Writing output");
    try
    {
        files.WriteSubClones(sample);
        files.WriteParentTree(lcaTree);
        files.WriteMullerDataFrames(aboveCutOff, connectedTree);
        files.WriteCCF(vaf, popSizes.Last().total);
        files.AddToSummary(resultSummary);
        files.StoreCopy(i);
        // files.WriteCopyNumbers(sample);
        // files.WriteRawData(random, sample, snps, simParams.IsFemale);
    }
    catch (Exception e)
    {
        Console.WriteLine($"Failed to write to disk with error: {e.Message}");
    }

    watch.Stop();
    Console.WriteLine($"Execution Time: {watch.ElapsedMilliseconds / 1000.0:F2}s\n");
}

files.CopySummary();
globalWatch.Stop();
Console.WriteLine($"Total time: {TimeSpan.FromMilliseconds(globalWatch.ElapsedMilliseconds)}");

return 0;