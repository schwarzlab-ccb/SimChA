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
    DivisionRate = 0.01f, 
    MutationRate = 0.01f, 
    IsFemale = true,
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
        [AbberationEnum.Chromothripsis] = 1f
    }
};

var simulator = new Simulator(simParams);
int stopCount = options.Value.StopCount;
int stepNo = 1;
while (simulator.Clones.Count < stopCount)
{
    Console.WriteLine($"Sim step {stepNo++:D3}, clones: {simulator.Clones.Count}");
    simulator.Step();
}

Console.WriteLine("Finished");
Console.WriteLine($"Total length is {ReferenceGenome.TotalLength(true)}");
Console.WriteLine($"Cell count {simulator.Clones.Sum(c => c.AliveCount)}");
int cutOffCount = simulator.Clones.Count(subClone => subClone.AliveCount >= options.Value.CutOff);
Console.WriteLine($"SubClone count {simulator.Clones.Count}. Above cutoff: {cutOffCount}");

try
{
    var files = new FileIO(options.Value.OutputPath);
    files.WriteSubClones(simulator.Clones, options.Value.CutOff);
    files.WriteParentGraph(simulator.Clones, options.Value.CutOff);
    files.WriteCopyNumbers(simulator.Clones, options.Value.CutOff);
} 
catch (Exception e) 
{
    Console.Write($"Failed to write to disk with error: {e.Message}");
}

var exampleSubClone = simulator.Clones.Where(subClone => subClone.AliveCount >= options.Value.CutOff).Last();
var rawData = new RawDataSingleSubclone(
    CopyNumbers.CalcCopyNumbers(exampleSubClone.Karyotype), 
    exampleSubClone.Karyotype.IsFemale, 
    1.0f, 
    2.0f);
rawData.CalcRawData();

string fullPathBAF = Path.Combine(Path.GetFullPath(options.Value.OutputPath), "SimCha_BAF.out");
Console.WriteLine($"Writing BAF for Subclone {exampleSubClone.CloneId} to file {fullPathBAF}");
using var outputFileBAF = new StreamWriter(fullPathBAF);
outputFileBAF.Write(rawData.BAFToTSV(true));

string fullPathlogR = Path.Combine(Path.GetFullPath(options.Value.OutputPath), "SimCha_logR.out");
Console.WriteLine($"Writing logR for Subclone {exampleSubClone.CloneId} to file {fullPathlogR}");
using var outputFilelogR = new StreamWriter(fullPathlogR);
outputFilelogR.Write(rawData.logRToTSV(true));

