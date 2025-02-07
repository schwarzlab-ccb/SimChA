using System.Diagnostics;
using System.Globalization;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Simulation;
using CommandLine;

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
var simParams = FileIO.ReadSimParams(options.ConfigFile);
var rnd = new Random(simParams.Seed);
var files = new FileIO(options.OutputPath);
var genRef = FileIO.GetGenRef(options.DataFolder, options.ShouldParseGenome);
var simulator = Factory.GetSimulator(rnd, genRef, simParams, selMode);
var samples = Factory.ReadSamples(rnd, genRef, simParams, options);

if (execMode != ExecMode.Profiles)
{
    Console.WriteLine("SIMULATION");
    foreach (var sample in samples)
    {
        simulator.SampleEvents(sample);
    }
}

Console.WriteLine("ANALYSIS");
files.WriteSimParams(simParams);
foreach (var sample in samples)
{
    int counter = 1;
    int total = sample.Clones.Count;
    foreach (var clone in sample.Clones)
    {
        Console.Write($"\rSample {sample.SampleId}. Clone {counter++}/{total}.".PadRight(80));
        sample.CloneStats[clone.CloneId] = CNProfile.GetCloneStats(sample, clone, genRef, simulator.FitnessParams, sample.Kars);
    }
}

Console.WriteLine("OUTPUT");
try
{
    files.WriteSamples(samples);
    files.WriteClones(samples);
    
    if (options.CalcConsistentCNs)
    {
        files.WriteConsistentCNs(genRef, samples);
    }
    if (ExecMode.Tree == execMode)
    {
        files.WriteTree(samples);
    }
    if (!options.LightweightOutput)
    {
        files.WriteCopyNumbers(genRef, samples);
        files.WriteKaryotypes(samples);
    }
    if (samples.Any(s => s.EventDescs.Any()))
    {
        files.WriteEvents(samples);
    }
    if (options.UseVariants)
    {
        files.WriteVCF(genRef, samples);
    }
    if (options.WriteFasta)
    {
        files.WriteFasta(genRef, samples);
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