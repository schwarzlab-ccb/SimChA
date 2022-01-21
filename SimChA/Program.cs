using SimChA.DataTypes;
using SimChA.Simulation;
using CommandLine;
using SimChA.Computation;
using SimChA.IO;

var options = Parser.Default.ParseArguments<CmdOptions>(args);
options.WithNotParsed(o =>
{
    Console.WriteLine("Exiting");
    o.ToList().ForEach(Console.Write);
    Environment.Exit(0);
});

var simParams = new SimParams
{
    DivisionRate = 0.02f,
    MutationRate = 0.01f,
    DeathRate = 0.01f,
    IsFemale = true,
    FitnessInc = 1f,
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

var simulator = new Simulator(simParams);
int stepNo = 1;
long pop = CellSampling.PopulationSize(simulator.Clones);
do
{
    Console.WriteLine($"Sim step {stepNo++:D3}, clones: {simulator.Clones.Count}, cells: {pop}");
    simulator.Step();
    pop = CellSampling.PopulationSize(simulator.Clones);
} while (pop < options.Value.StopCount && pop > 0);

Console.WriteLine("Finished");
Console.WriteLine($"Total length is {ReferenceGenome.TotalLength(true)}");
Console.WriteLine($"Cell count {simulator.Clones.Sum(c => c.AliveCount)}");
// var sample = CellSampling.SampleCells(simulator.Clones, 10000);
float cutOff = pop * options.Value.CutOff;
var sample = simulator.Clones.Where(sc => sc.AliveCount >= cutOff).ToList();
Console.WriteLine($"SubClone count {simulator.Clones.Count}. Above cutoff: { sample.Count }");
// snps are shared between all subclones and therefore are created only once
var snps = SNPs.CreateSNPs(simParams.IsFemale, 100);
var parentMap = TreeBuilder.CreateParentMap(simulator.Clones);
var parentTree = TreeBuilder.BuildTree(parentMap, sample);
try
{
    var files = new FileIO(options.Value.OutputPath);
    files.WriteSubClones(sample);
    files.WriteParentTree(parentTree);
    files.WriteCopyNumbers(sample);
    files.WriteRawData(sample, snps, simParams.IsFemale);
    files.WriteMullerDataFrames(sample);
} 
catch (Exception e) 
{
    Console.WriteLine($"Failed to write to disk with error: {e.Message}");
}
