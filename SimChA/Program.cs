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
    Environment.Exit(0);
});

var watch = new System.Diagnostics.Stopwatch();
watch.Start();

int seed = options.Value.Seed >= 0 ? options.Value.Seed : new Random().Next();
var simParams = new SimParams
{
    Seed = seed,
    PopLimit = options.Value.StopCount,
    CutOff = options.Value.CutOff,
    IsFemale = true,
    DivisionRate = 0.02f,
    MutationRate = 0.01f,
    DriverProb = .01f,
    FitnessIncMu = .01f,
    FitnessIncSigma = 1, // Multiplication of the Division rate
    DeathRate = 0.99f, // Multiplication of the Division rate
    SplitRate = 0.0f,
    Confinement = 0.0f,
    DecayRate = 0.0f,
    InitialPop = 25,
    
    StepLimit = 10000,
    IsMultiplicative = false,
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

var random = new Random(simParams.Seed);
var simulator = new Simulator(simParams, random);
int stepNo = 0;
var popSizes = new List<(long total, long alive)>
{
    (CellSampling.PopulationSize(simulator.Populations), CellSampling.AliveCount(simulator.Populations))
};
Console.WriteLine($"Sim with seed {seed}, genome length  {ReferenceGenome.TotalLength(true)}");
int clones;
do
{
    clones = simulator.FlatPops.Count();
    Console.WriteLine($"Step: {++stepNo:D3}, " +
                      $"populations: {simulator.Populations.Count}, " +
                      $"subClones: {clones}, " +
                      $"cells: { popSizes.Last().total }, " +
                      $"alive: { popSizes.Last().alive }");
    simulator.Step();
    popSizes.Add((CellSampling.PopulationSize(simulator.Populations), CellSampling.AliveCount(simulator.Populations)));
} while (clones < simParams.PopLimit && popSizes.Last().alive > 0 && stepNo < simParams.StepLimit);


var cutOff = popSizes.Select(pair => (long) Math.Ceiling(pair.alive * simParams.CutOff)).ToList();
var aboveCutOff = simulator.FlatPops.Where(sc 
    => Enumerable.Range(0, popSizes.Count).Any(g => cutOff[g] <= sc.AliveAtGen(g))).ToList();
var lcaTree = LCATreeBuilder.Builtree(simulator.FlatPops, aboveCutOff);
var connectedTree = ConnectedTreeBuilder.BuildTree(simulator.FlatPops, aboveCutOff);
var treeNodes = lcaTree.Nodes.Select(n => n.Id).ToList();
var sample = simulator.FlatPops.Where(sc => treeNodes.Contains(sc.CloneId)).ToList();
var snps = SNPBuilder.CreateSNPs(random, simParams.IsFemale, 100); // snps are shared between all subclones and therefore are created only once
var vaf = TreeAnalysis.ComputeVAF(connectedTree);
(int nodeCount, int leafCount, int depth, float branching) = TreeAnalysis.ComputeTreeSize(connectedTree);
Console.WriteLine($"SubClone count {simulator.FlatPops.Count()}. Above cutoff: { sample.Count }");
Console.WriteLine($"nodeCount: {nodeCount}, leafCount: {leafCount}, depth: {depth}, branching: {branching:F2} ");

try
{
    var files = new FileIO(options.Value.OutputPath);
    files.WriteSimParams(simParams);
    files.WriteSubClones(sample);
    files.WriteParentTree(lcaTree);
    files.WriteMullerDataFrames(aboveCutOff, connectedTree);
    files.WriteCopyNumbers(sample);
    files.WriteVAF(vaf, popSizes.Last().total);
    files.WriteRawData(random, sample, snps, simParams.IsFemale);
} 
catch (Exception e) 
{
    Console.WriteLine($"Failed to write to disk with error: {e.Message}");
}

watch.Stop();
Console.WriteLine($"Execution Time: {(watch.ElapsedMilliseconds/1000.0):F2}s");