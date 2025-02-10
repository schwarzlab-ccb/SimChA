using System.Diagnostics;
using System.Globalization;
using SimChA.Computation;
using SimChA.IO;
using SimChA.Simulation;
using CommandLine;
using SimChA.Data;

// Configuration
Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
var cmdOptions = Parser.Default.ParseArguments<CmdOptions>(args);
cmdOptions.WithNotParsed(o =>
{
    Console.WriteLine("Exiting");
    o.ToList().ForEach(Console.WriteLine); // Write out errors
    Environment.Exit(1);
});
var watch = new Stopwatch();
watch.Start();

Console.WriteLine("INPUT");
var options = cmdOptions.Value;
var execMode = options.ExecMode;
var selMode = options.SelectionMode;
var config = FileIO.ReadSimParams(options.ConfigFile);
var rnd = new Random(config.Seed);
var files = new FileIO(options.OutputPath);
var genRef = FileIO.GetGenRef(options.DataFolder, options.ShouldParseGenome);
var samples = Factory.ReadSamples(rnd, genRef, config, options);
var simulator = Factory.GetSimulator(rnd, genRef, config, selMode);
var clones = new Dictionary<string, List<Clone>>();
    
if (options.Simulate)
{
    Console.WriteLine("SIMULATION");
    foreach (var sample in samples)
    {
        clones[sample.SampleId] = simulator.Simulate(sample);
    }
}

Console.WriteLine("ANALYSIS");
var cloneList = new List<CloneStat>();
        
foreach ((string sampleId, var subClones) in clones)
{
    int counter = 1;
    int total = clones.Values.Count;
    foreach (var clone in subClones)
    {
        Console.Write($"\rSample {sampleId}. Clone {counter++}/{total}.".PadRight(80));
        var sampleStats = subClones.Select(s 
            => CNProfile.GetCloneStats(sampleId, clone, genRef, config.FitParams));
        cloneList.AddRange(sampleStats);
    }
}

Console.WriteLine("OUTPUT");
files.WriteSimParams(config);
try
{
    files.WriteSamples(samples);
    files.WriteClones(cloneList);
    
    if (options.CalcConsistentCNs)
    {
        files.WriteConsistentCNs(genRef, clones);
    }
    if (ExecMode.Tree == execMode)
    {
        files.WriteTree(samples);
    }
    if (!options.LightweightOutput)
    {
        files.WriteCopyNumbers(genRef, clones);
        files.WriteKaryotypes(clones);
    }
    if (options.Simulate)
    {
        files.WriteEvents(clones);
    }
    if (options.UseVariants)
    {
        files.WriteVCF(genRef, clones);
    }
    if (options.WriteFasta)
    {
        files.WriteFasta(genRef, clones);
    }
}
catch (Exception e)
{
    Console.WriteLine($"Failed with exception {e.Message}. Stack: {e.StackTrace}");
    return e.HResult;
}

watch.Stop();
Console.WriteLine($"Total time: {TimeSpan.FromMilliseconds(watch.ElapsedMilliseconds)}");

return 0;