using SimChA.DataTypes;
using SimChA.Simulation;
using CommandLine;
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

string fullPathSubclones = Path.Combine(Path.GetFullPath(options.Value.OutputPath), "SimCha_subclones.out");
Console.WriteLine($"Writing subclones to file {fullPathSubclones}");
using var outputFileSubclones = new StreamWriter(fullPathSubclones);
foreach (var subClone in simulator.Clones.Where(subClone => subClone.AliveCount >= options.Value.CutOff))
{
    outputFileSubclones.Write(subClone);
}

string fullPathCopyNumbers = Path.Combine(Path.GetFullPath(options.Value.OutputPath), "SimCha_copynumbers.out");
Console.WriteLine($"Writing CopyNumbers to file {fullPathCopyNumbers}");
using var outputFileCopyNumbers = new StreamWriter(fullPathCopyNumbers);
outputFileCopyNumbers.Write("sample_id\tchrom\tstart\tend\tcn_a\tcn_b\n");
foreach (var subClone in simulator.Clones.Where(subClone => subClone.AliveCount >= options.Value.CutOff))
{
    var copynumbers = new CopyNumbers(subClone.Karyotype);
    outputFileCopyNumbers.Write(copynumbers.ToTSV(subClone.CloneId.ToString(), false));
}

var exampleSubClone = simulator.Clones.Where(subClone => subClone.AliveCount >= options.Value.CutOff).Last();
var rawData = new RawDataSingleSubclone(new CopyNumbers(exampleSubClone.Karyotype), 1.0f, 2.0f);

string fullPathBAF = Path.Combine(Path.GetFullPath(options.Value.OutputPath), "SimCha_BAF.out");
Console.WriteLine($"Writing BAF for Subclone {exampleSubClone.CloneId} to file {fullPathBAF}");
using var outputFileBAF = new StreamWriter(fullPathBAF);
outputFileBAF.Write(rawData.BAFToTSV(true));

string fullPathlogR = Path.Combine(Path.GetFullPath(options.Value.OutputPath), "SimCha_logR.out");
Console.WriteLine($"Writing logR for Subclone {exampleSubClone.CloneId} to file {fullPathlogR}");
using var outputFilelogR = new StreamWriter(fullPathlogR);
outputFilelogR.Write(rawData.logRToTSV(true));

