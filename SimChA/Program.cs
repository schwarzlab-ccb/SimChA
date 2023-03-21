using System.Diagnostics;
using CommandLine;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Simulation;

var watch = new Stopwatch();
watch.Start();

var options = Parser.Default.ParseArguments<CmdOptions>(args);
options.WithNotParsed(o =>
{
    Console.WriteLine("Exiting");
    o.ToList().ForEach(Console.WriteLine); // Write out errors
    Environment.Exit(1);
});

// MCMC simulation vs Adam and Tom's CNP stuff:
bool mcmc_on = options.Value.MCMC_ON;

if (mcmc_on)
{
    Console.WriteLine("Hello, world :)");
}
else
{
    SimParams simParams;
    string configFile = options.Value.ConfigFile;
    if (configFile != "")
    {
        simParams = FileIO.ReadSimParams(configFile);
    }
    else
    {
        int seed = new Random().Next();
        var fitness = new FitnessParams(1, 1, 1);
        simParams = new SimParams(seed, true, 1, Distribution.Uniform, GenomeAssembly.hg38, fitness, null);
    }

    HGRef.Assembly = simParams.Assembly;
    var rnd = new Random(simParams.Seed);
    var files = new FileIO(options.Value.OutputPath);
    var geneLists = FileIO.ReadGeneLists(options.Value.GenesFolder, simParams.SexXX, HGRef.Assembly);
    files.WriteSimParams(simParams);

    // Obtain clones
    var newick = "";
    List<Clone> clones = new();
    List<CNEvent> cnEvents = new();
    if (options.Value.CNProfiles != "")
    {
        var profiles = FileIO.ReadProfiles(options.Value.CNProfiles);
        clones = Simulator.ClonesFromProfiles(profiles);
    }
    else
    {
        Parsers.ValidateSignatures(simParams.Signatures);
        Console.WriteLine("Computing mutations.");
        var simulator = new Simulator(rnd, simParams, geneLists);

        if (options.Value.NewickFile != "")
        {
            newick = FileIO.ReadNewick(options.Value.NewickFile);
            clones = Parsers.ParseNewick(newick, simParams.SexXX);
        }
        else
        {
            clones = Simulator.MakeClones(rnd, options.Value.Repeats, simParams.SexXX, simParams.EventCount, simParams.Distribution);   
        }
        cnEvents = simulator.ApplyEvents(clones[0], clones);
        clones = clones.Where(c => c.CloneId != 0).ToList();
    }

    // Fitness data
    Console.WriteLine("Computing fitness.");
    var results = new List<ProfileStats>();
    var counter = 0;
    foreach (var clone in clones)
    {
        Console.Write($"\r{++counter}/{clones.Count} samples processed.");
        var profileStats = CNProfile.GetProfileStats(clone, geneLists, simParams.Fitness);
        results.Add(profileStats);
    }

    Console.WriteLine();
    try
    {
        files.WriteClones(clones);
        files.WriteCopyNumbers(clones, simParams.SexXX);
        files.WriteSampleFitness(results);
        if (cnEvents.Any())
        {
            files.WriteEvents(cnEvents);
        }
        if (newick != "")
        {
            files.WriteNewick(newick);
        }
    }
    catch (Exception e)
    {
        Console.WriteLine($"Failed with exception {e.Message}. Stack: {e.StackTrace}");
        return e.HResult;
    }
}


watch.Stop();
Console.WriteLine($"Total time: {TimeSpan.FromMilliseconds(watch.ElapsedMilliseconds)}");

return 0;