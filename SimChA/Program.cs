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
    Environment.Exit(1);
});

SimParams simParams;
if (options.Value.ConfigFile != "")
{
    simParams = FileIO.SimParamsFromFile(options.Value.ConfigFile);
}
else
{
    simParams = new SimParams
    {
        Checkpoints = false,
        // Function
        FitnessAcc = FitnessAccType.Mul,
        FitnessDist = FitnessSampleType.Constant,
        Seed = new Random().Next(),
        // Experiment
        PopLimit = 1_03_741_824,
        StepLimit = 1_000_000,
        CutOff = 0.01f,
        Repeats = 1,
        InitPop = 100,

        // Model
        Turnover = 0.01,
        
        MutationProb = 0.0001,
        FitnessMean = .1,
        Confinement = .1,
        
        InitMut = 1,
    };
}

var random = new Random(simParams.Seed);
FileIO files;
bool isRepeated = simParams.Repeats > 1;
try
{
    files = new FileIO(options.Value.OutputPath, isRepeated);
    files.WriteSimParams(simParams);
}
catch (Exception e)
{
    Console.WriteLine($"Failed to write to disk with error: {e.Message}");
    return 2;
}

try
{
    var globalWatch = new System.Diagnostics.Stopwatch();
    globalWatch.Start();

    for (int repeatId = 0; repeatId < simParams.Repeats; repeatId++)
    {
        var watch = new System.Diagnostics.Stopwatch();
        watch.Start();

        // Simulation
        int checkpointId = 0;
        var simulator = new Simulator(simParams, random);
        var checkpoints = Utility.CreateCheckpoints(simParams);
        var popSizes = new List<(long total, long alive, long lost)> {CellSampling.PopState(simulator.Clones)};
        var EndCondFunc = () =>
            !(popSizes.Last().total <= simParams.PopLimit
              && popSizes.Last().alive > 0
              && simulator.StepNo < simParams.StepLimit);
        do
        {
            simulator.Step();
            popSizes.Add(CellSampling.PopState(simulator.Clones));
            if (popSizes.Last().alive <= 0 && popSizes.Last().total < simParams.InitPop)
            {
                break;
            }

            Console.Write(($"Sim: {repeatId + 1}/{simParams.Repeats}, " +
                           $"step: {simulator.StepNo:D3}, " +
                           $"subClones: {simulator.Clones.Count}, " +
                           $"alive SC: {simulator.AliveSC}, " +
                           $"cells: {popSizes.Last().total:N0}, " +
                           $"alive: {popSizes.Last().alive:N0}, " +
                           $"lost: {popSizes.Last().lost:N0}").PadRight(160) +
                          (options.Value.Newline ? "\n" : "\r"));

            if (EndCondFunc() || checkpoints.Any() && popSizes.Last().total > checkpoints.First())
            {
                // Analysis
                var cutOff = popSizes.Select(pair => (long) Math.Ceiling(pair.alive * simParams.CutOff)).ToList();
                var aboveCutOff = simulator.Clones
                    .Where(sc => Enumerable.Range(0, popSizes.Count).Any(g => cutOff[g] <= sc.AliveAtGen(g))).ToList();
                var lcaTree = LCATreeBuilder.Builtree(simulator.Clones, aboveCutOff);
                var connectedTree = ConnectedTreeBuilder.BuildTree(simulator.Clones, aboveCutOff);
                var treeNodes = lcaTree.Nodes.Select(n => n.Id).ToList();
                var sample = simulator.Clones.Where(sc => treeNodes.Contains(sc.CloneId)).ToList();

                // Summary
                if (checkpoints.Any())
                {
                    checkpointId++;
                    checkpoints.RemoveAt(0);
                }

                string time = TimeSpan.FromMilliseconds(watch.ElapsedMilliseconds).ToString();
                var result = new ResultSummary(repeatId, checkpointId, connectedTree, aboveCutOff, 
                    simulator.Clones.Count, simulator.AliveSC, sample.Count, simulator.StepNo, popSizes,
                    time);
                files.AddToSummary(result);

                // Result
                if (EndCondFunc())
                {
                    var vaf = TreeAnalysis.ComputeVAF(connectedTree);
                    files.WriteFinalOutput(repeatId, sample, lcaTree, aboveCutOff, connectedTree, vaf,
                        popSizes.Last().total);
                    Console.WriteLine($"Sim: {repeatId + 1}/{simParams.Repeats} result:".PadRight(160));
                    Console.WriteLine(result.ToText());
                }
            }
        } while (!EndCondFunc());

        watch.Stop();
        Console.WriteLine($"Execution Time: {TimeSpan.FromMilliseconds(watch.ElapsedMilliseconds)}");
        Console.WriteLine(string.Join("", Enumerable.Repeat("*", 100)));
    }

    files.CopySummary();
    globalWatch.Stop();

    Console.WriteLine($"Total time: {TimeSpan.FromMilliseconds(globalWatch.ElapsedMilliseconds)}");
}
catch (Exception e)
{
    Console.WriteLine($"Failed with exception {e.Message}");
}

return 0;