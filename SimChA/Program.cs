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

int seed = options.Value.Seed >= 0 ? options.Value.Seed : new Random().Next();
var simParams = new SimParams
{
    Seed = seed,
    PopLimit = options.Value.StopCount,
    CutOff = options.Value.CutOff,
    IsFemale = true,
    DivisionRate = 0.001f,
    MutationRate = 0.00f,
    DriverProb = 0.0f,
    FitnessIncMu = .01f,
    FitnessIncSigma = 1, // Multiplication of the Division rate
    DeathRate = 1f, // Multiplication of the Division rate
    SplitRate = 0.0f,
    Confinement = 0.0f,
    DecayRate = 0.0f,
    InitialPop = 100,
    
    StepLimit = 10_000,
    IsMultiplicative = false,
    // AberrationRates =
    // {
    //     [AberrationEnum.InternalDeletion] = 50f,
    //     [AberrationEnum.InternalDuplication] = 50f,
    //     [AberrationEnum.Translocation] = 20f,
    //     [AberrationEnum.TailDeletion] = 15f,
    //     [AberrationEnum.BreakageFusionBridge] = 10f,
    //     [AberrationEnum.Inversion] = 10f,
    //     [AberrationEnum.Missegregation] = 5f,
    //     [AberrationEnum.Duplication] = 5f,
    //     [AberrationEnum.Chromothripsis] = 1f,
    //     [AberrationEnum.WholeGenomeDoubling] = 1f
    // }
};

var random = new Random(simParams.Seed);
FileIO files;
try
{
    files = new FileIO(options.Value.OutputPath);
}
catch (Exception e)
{
    Console.WriteLine($"Failed to write to disk with error: {e.Message}");
    return 2;
}

var globalWatch = new System.Diagnostics.Stopwatch();
globalWatch.Start();

for (int i = 0; i < options.Value.Reps; i++) {
    var watch = new System.Diagnostics.Stopwatch();
    watch.Start();
    
    Console.WriteLine(string.Join("", Enumerable.Repeat("*" , 100 )));
    Console.WriteLine($"* Simulation {i}/{options.Value.Reps}");
    Console.WriteLine($"* Sim with seed {seed}, genome length  {ReferenceGenome.TotalLength(true)}");
    
    // Simulation
    var simulator = new Simulator(simParams, random);
    int stepNo = 0;
    var popSizes = new List<(long total, long alive)>();
    popSizes.Add((CellSampling.PopulationSize(simulator.Populations), CellSampling.AliveCount(simulator.Populations)));
    int clones;
    do
    {
        clones = simulator.FlatPops.Count();
        Console.Write($"Step: {++stepNo:D3}, " +
                      $"populations: {simulator.Populations.Count}, " +
                      $"subClones: {clones}, " +
                      $"cells: {popSizes.Last().total}, " +
                      $"alive: {popSizes.Last().alive}" +
                      "\r");
        simulator.Step();
        popSizes.Add((CellSampling.PopulationSize(simulator.Populations), CellSampling.AliveCount(simulator.Populations)));
    } while (clones < simParams.PopLimit && popSizes.Last().alive > 0 && stepNo < simParams.StepLimit);
    
    // Analysis
    var cutOff = popSizes.Select(pair => (long)Math.Ceiling(pair.alive * simParams.CutOff)).ToList();
    var aboveCutOff = simulator.FlatPops.Where(sc
        => Enumerable.Range(0, popSizes.Count).Any(g => cutOff[g] <= sc.AliveAtGen(g))).ToList();
    var lcaTree = LCATreeBuilder.Builtree(simulator.FlatPops, aboveCutOff);
    var connectedTree = ConnectedTreeBuilder.BuildTree(simulator.FlatPops, aboveCutOff);
    var treeNodes = lcaTree.Nodes.Select(n => n.Id).ToList();
    var sample = simulator.FlatPops.Where(sc => treeNodes.Contains(sc.CloneId)).ToList();
    var vaf = TreeAnalysis.ComputeVAF(connectedTree);
    // var snps = SNPBuilder.CreateSNPs(random, simParams.IsFemale, 100); // snps are shared between all subclones and therefore are created only once
    
    // Summary
    ResultSummary resultSummary = new ();
    (resultSummary.nodeCount, resultSummary.leafCount, resultSummary.treeDepth, resultSummary.branching) 
        = TreeAnalysis.ComputeTreeSize(connectedTree);
    resultSummary.subcloneTotal = clones;
    resultSummary.subcloneSelect = sample.Count;
    resultSummary.generations = stepNo;
    resultSummary.aliveCount = popSizes.Last().alive;
    resultSummary.totalCount = popSizes.Last().total;
    Console.WriteLine($"SubClone count {resultSummary.subcloneTotal}. Above cutoff: {resultSummary.subcloneSelect}");
    Console.WriteLine(ResultSummary.Header());
    Console.WriteLine(resultSummary.ToString());

    // Output
    if (i == 0)
    {
        try
        {
            files.WriteSimParams(simParams);
            files.WriteSubClones(sample);
            files.WriteParentTree(lcaTree);
            files.WriteMullerDataFrames(aboveCutOff, connectedTree);
            files.WriteCCF(vaf, popSizes.Last().total);
            files.CreateSummary();
            // files.WriteCopyNumbers(sample);
            // files.WriteRawData(random, sample, snps, simParams.IsFemale);
        }
        catch (Exception e)
        {
            Console.WriteLine($"Failed to write to disk with error: {e.Message}");
        }
    }
    files.AddToSummary(resultSummary);
    
    
    watch.Stop();
    Console.WriteLine($"Execution Time: {watch.ElapsedMilliseconds/1000.0:F2}s\n");
}

globalWatch.Stop();
Console.WriteLine($"Total time: {TimeSpan.FromMilliseconds(globalWatch.ElapsedMilliseconds)}");

return 0;