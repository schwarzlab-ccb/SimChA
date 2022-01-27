using System.Collections;
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

var simParams = new SimParams
{
    IsFemale = true,
    DivisionRate = 0.025f,
    MutationRate = 0.01f,
    DriverProb = 0.1f,
    DeathRate = 0.01f,
    SplitRate = 0.0005f,
    DivisionSlowDown = 0.005,
    FitnessInc = 1.2f,
    InitialPop = 100,
    AbberationRates =
    {
        [AbberationEnum.InternalDeletion] = 50f,
        [AbberationEnum.InternalDuplication] = 50f,
        [AbberationEnum.Translocation] = 20f,
        [AbberationEnum.TailDeletion] = 15f,
        [AbberationEnum.BreakageFusionBridge] = 10f,
        [AbberationEnum.Inversion] = 10f,
        [AbberationEnum.Missegregation] = 5f,
        [AbberationEnum.Duplication] = 5f,
        [AbberationEnum.Chromothripsis] = 1f,
        [AbberationEnum.WholeGenomeDoubling] = 1f
    }
};

int seed = options.Value.Seed >= 0 ? options.Value.Seed : new Random().Next();
var random = new Random(seed);
var simulator = new Simulator(simParams, random);
int stepNo = 0;
var popSizes = new List<long> { CellSampling.PopulationSize(simulator.Populations) };
do
{
    Console.WriteLine($"Sim step {++stepNo:D3}, " +
                      $"populations: {simulator.Populations.Count}, " +
                      $"subclones: {simulator.FlatPops.Count()}, " +
                      $"cells: { popSizes[^1] }");
    simulator.Step();
    popSizes.Add(CellSampling.PopulationSize(simulator.Populations));
} while (popSizes[^1] < options.Value.StopCount && popSizes[^1] > 0);

Console.WriteLine("Finished");
Console.WriteLine($"Seed used was {seed}");
Console.WriteLine($"Total length is {ReferenceGenome.TotalLength(true)}");
Console.WriteLine($"Cell count {popSizes[^1]}");

// snps are shared between all subclones and therefore are created only once
var snps = SNPBuilder.CreateSNPs(random, simParams.IsFemale, 100);
var cutOff = popSizes.Select(l => (long) Math.Ceiling(l * options.Value.CutOff)).ToList();
var aboveCutOff = simulator.FlatPops.Where(sc 
    => Enumerable.Range(sc.FirstGen, popSizes.Count - sc.FirstGen).Any(g => cutOff[g] <= sc.PopAtGeneration(g))
    ).ToList();
var lcaTree = LCATreeBuilder.Builtree(simulator.FlatPops, aboveCutOff);
var connectedTree = ConnectedTreeBuilder.BuildTree(simulator.FlatPops, aboveCutOff);
var treeNodes = lcaTree.Nodes.Select(n => n.Id).ToList();
var sample = simulator.FlatPops.Where(sc => treeNodes.Contains(sc.CloneId)).ToList();
Console.WriteLine($"SubClone count {simulator.Populations.Count}. Above cutoff: { sample.Count }");

try
{
    var files = new FileIO(options.Value.OutputPath);
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
