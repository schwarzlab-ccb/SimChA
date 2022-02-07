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

int seed = options.Value.Seed >= 0 ? options.Value.Seed : new Random().Next();
var simParams = new SimParams
{
    Seed = seed,
    PopLimit = options.Value.StopCount,
    CutOff = options.Value.CutOff,
    IsFemale = true,
    DivisionRate = 0.025f,
    MutationRate = 0.01f,
    DriverProb = 0.1f,
    DeathRate = 0.01f,
    SplitRate = 0.0f,
    DivisionSlowDown = 0.001f,
    DecayRate = 0.01f,
    FitnessIncMu = 1.05f,
    FitnessIncSigma = .1f,
    InitialPop = 1000,
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

var random = new Random(seed);
var simulator = new Simulator(simParams, random);
int stepNo = 0;
var popSizes = new List<(long, long)> { (CellSampling.PopulationSize(simulator.Populations), CellSampling.AliveCount(simulator.Populations)) };
Console.WriteLine($"Sim with seed {seed}, genome length  {ReferenceGenome.TotalLength(true)}");
do
{
    Console.WriteLine($"Step: {++stepNo:D3}, " +
                      $"populations: {simulator.Populations.Count}, " +
                      $"subClones: {simulator.FlatPops.Count()}, " +
                      $"cells: { popSizes.Last().Item1 }, " +
                      $"alive: { popSizes.Last().Item2 }");
    simulator.Step();
    popSizes.Add((CellSampling.PopulationSize(simulator.Populations), CellSampling.AliveCount(simulator.Populations)));
} while (popSizes.Last().Item1 < simParams.PopLimit && popSizes.Last().Item2 > 0);



var cutOff = popSizes.Select(pair => (long) Math.Ceiling(pair.Item2 * simParams.CutOff)).ToList();
var aboveCutOff = simulator.FlatPops.Where(sc 
    => Enumerable.Range(0, popSizes.Count).Any(g => cutOff[g] <= sc.AliveAtGen(g))).ToList();
var lcaTree = LCATreeBuilder.Builtree(simulator.FlatPops, aboveCutOff);
var connectedTree = ConnectedTreeBuilder.BuildTree(simulator.FlatPops, aboveCutOff);
var treeNodes = lcaTree.Nodes.Select(n => n.Id).ToList();
var sample = simulator.FlatPops.Where(sc => treeNodes.Contains(sc.CloneId)).ToList();
var snps = SNPBuilder.CreateSNPs(random, simParams.IsFemale, 100); // snps are shared between all subclones and therefore are created only once
Console.WriteLine($"SubClone count {simulator.FlatPops.Count()}. Above cutoff: { sample.Count }");

try
{
    var files = new FileIO(options.Value.OutputPath);
    files.WriteSimParams(simParams);
    files.WriteSubClones(sample);
    files.WriteParentTree(lcaTree);
    files.WriteMullerDataFrames(aboveCutOff, connectedTree);
    files.WriteCopyNumbers(sample);
    files.WriteRawData(random, sample, snps, simParams.IsFemale);
} 
catch (Exception e) 
{
    Console.WriteLine($"Failed to write to disk with error: {e.Message}");
}
