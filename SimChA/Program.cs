using System.Diagnostics;
using System.Globalization;
using SimChA.IO;
using SimChA.Simulation;
using SimChA.Data;
using CommandLine;

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

// Input
var options = cmdOptions.Value;
var selMode = options.SelectionMode;
var config = FileIO.ReadSimChAConfig(options.ConfigFile);
var rnd = new Random(config.ChAParams.Seed);
var files = new FileIO(options.OutputPath);
var genRef = FileIO.GetGenRef(options.DataFolder, options.ShouldParseGenome);
var simulator = Factory.GetSimulator(rnd, genRef, config, selMode);

// Construct samples
var samples = new List<Sample>();
if (options.Simulate)
{
    var validSigs = Factory.ValidateSignatures(config.Signatures);
    if (options.ExecMode == ExecMode.Repeats)
    {
        Console.WriteLine($"Creating {options.Repeats} samples:");
        for (int i = 0; i < options.Repeats; i++)
        {
            string sampleId = $"Sample_{i}";
            var node = new CTreeNode(sampleId, sampleId, -1, -1);
            var tree = new List<CTreeNode> {node};
            var newSample = simulator.Simulate(node, tree, validSigs);
            samples.Add(newSample[0]);
        }
    }
    else
    {
        var (root, tree) = FileIO.ReadCloneTree(options.CloneTreeFile, options.MHMode);
        Console.WriteLine($"Creating {tree.Count} samples from a tree:");
        samples = simulator.Simulate(root, tree, validSigs);
    }
}
else
{
    Console.WriteLine("Reading samples from data:");
    samples = FileIO.ReadProfiles(genRef, options.CNProfiles, config.ChAParams.AutosomesOnly);
}

// Score samples
var sampleStats = new List<SampleStat>();
for (int i = 0; i < samples.Count; i++)
{
    Console.Write($"\rAnalyzing sample {i + 1}/{samples.Count}.".PadRight(80));
    var sampleStat = SampleStat.GetSampleStat(samples[i], genRef, config.FitParams);
    sampleStats.Add(sampleStat);
}

// Output
Console.WriteLine();
files.WriteSimParams(config);
try
{
    files.WriteSamples(sampleStats);
    
    if (options.CalcConsistentCNs)
    {
        files.WriteConsistentCNs(genRef, samples);
    }
    if (!options.LightweightOutput)
    {
        files.WriteCopyNumbers(genRef, samples);
        files.WriteKaryotypes(samples);
    }
    if (options.Simulate)
    {
        files.WriteEvents(samples);
    }
    if (options.UseVariants)
    {
        files.WriteVCF(genRef, samples);
    }
    if (options.WriteFasta)
    {
        files.WriteFasta(samples);
    }
}
catch (Exception e)
{
    Console.WriteLine($"Failed to write data with exception {e.Message}. Stack: {e.StackTrace}");
    return e.HResult;
}

watch.Stop();
Console.WriteLine($"Total time: {TimeSpan.FromMilliseconds(watch.ElapsedMilliseconds)}");
return 0;