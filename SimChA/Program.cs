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

// TestCopyNumbers

string fullPathCopyNumbers = Path.Combine(Path.GetFullPath(options.Value.OutputPath), "SimCha_copynumbers.out");
Console.WriteLine($"Writing CopyNumbers to file {fullPathCopyNumbers}");
using var outputFileCopyNumbers = new StreamWriter(fullPathCopyNumbers);
outputFileCopyNumbers.Write("sample_id\tchrom\tstart\tend\tcn_a\tcn_b\n");
foreach (var subClone in simulator.Clones.Where(subClone => subClone.AliveCount >= options.Value.CutOff))
{
    var copynumbers = new CopyNumbers(subClone.Karyotype);
    outputFileCopyNumbers.Write(copynumbers.ToTSV(subClone.CloneId.ToString(), false));
}
