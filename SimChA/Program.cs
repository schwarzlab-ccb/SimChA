using System.Diagnostics;
using CommandLine;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Simulation;

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
        Checkpoints = true,
        // Function
        FitnessAcc = FitnessAccType.Mul,
        FitnessDist = FitnessSampleType.Constant,
        FitnessEffect = FitnessEffectType.Death,
        Seed = new Random().Next(),
        // Experiment
        MinPop = 1000,
        MaxPop = 1_08_576_000,
        MaxSteps = 1_000_000,
        CutOff = 0.001f,
        Repeats = 1,

        // Model
        Turnover = 0.01,
        MutationProb = 0.0001,

        FitnessMean = .125,
        Confinement = .1,

        // Initialization
        StartMut = 1,
        StartPop = 1
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
    var globalWatch = new Stopwatch();
    globalWatch.Start();

    int tryNo = 0;
    for (int repeatId = 0; repeatId < simParams.Repeats; repeatId++)
    {
        var watch = new Stopwatch();
        watch.Start();

        // Simulation
        int checkpointId = 0;
        var simulator = new Simulator(simParams, random);
        var checkpoints = Utility.CreateCheckpoints(simParams);
        var popSizes = new List<(long total, long alive, long necro)> { CellSampling.PopState(simulator.Clones) };
        var EndCondFunc = () =>
            !(popSizes.Last().total <= simParams.MaxPop
              && popSizes.Last().alive > 0
              && simulator.StepNo < simParams.MaxSteps);
        do
        {
            simulator.Step();
            popSizes.Add(CellSampling.PopState(simulator.Clones));
            Console.Write(($"sim: {repeatId + 1}.{tryNo}/{simParams.Repeats}, " +
                           $"step: {simulator.StepNo:D3}, " +
                           $"SC_total: {simulator.Clones.Count}, " +
                           $"SC_alive: {simulator.AliveSC}, " +
                           $"C_total: {popSizes.Last().total:N0}, " +
                           $"C_alive: {popSizes.Last().alive:N0}, " +
                           $"C_necro: {popSizes.Last().necro:N0}").PadRight(160) +
                          (options.Value.Newline ? "\n" : "\r"));

            if ((EndCondFunc() && popSizes.Last().total >= simParams.MinPop)
                || (checkpoints.Any() && popSizes.Last().total > checkpoints.First()))
            {
                // Analysis
                double cutOff = popSizes.Last().alive * simParams.CutOff;
                var aboveCutOff = simulator.Clones.Where(sc => sc.AliveCount > cutOff).ToList();
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
                var result = new ResultSummary(repeatId, checkpointId, simulator.StepNo, time,
                    connectedTree, aboveCutOff, simulator.Clones, popSizes);
                files.AddToSummary(result);

                // Result
                if (EndCondFunc() && popSizes.Last().total >= simParams.MinPop)
                {
                    var vaf = TreeAnalysis.ComputeVAF(connectedTree);
                    files.WriteSubClones(sample);
                    files.WriteParentTree(lcaTree);
                    files.WriteCCF(vaf, popSizes.Last().total);
                    var cutOffList = popSizes.Select(pair => pair.alive * 0.01).ToList();
                    int firstPop = cutOffList.FindIndex(minPop => minPop > 0);
                    var mullerPops = simulator.Clones.Where(sc =>
                        sc.FirstGen <= firstPop || Enumerable.Range(firstPop, popSizes.Count)
                            .Any(g => cutOffList[g] <= sc.AliveAtGen(g))).ToList();
                    var mullerTree = ConnectedTreeBuilder.BuildTree(simulator.Clones, mullerPops);
                    files.WriteMullerDataFrames(mullerPops, mullerTree);
                    files.StoreCopy(repeatId);
                    Console.WriteLine($"Sim: {repeatId + 1}.{tryNo}/{simParams.Repeats} result:".PadRight(160));
                    Console.WriteLine(result.ToText());
                }
            }
        } while (!EndCondFunc());

        // Skip on failure
        if (popSizes.Last().total < simParams.MinPop)
        {
            tryNo++;
            repeatId--;
        }
        else
        {
            tryNo = 0;
            watch.Stop();
            Console.WriteLine($"Execution Time: {TimeSpan.FromMilliseconds(watch.ElapsedMilliseconds)}");
            Console.WriteLine(string.Join("", Enumerable.Repeat("*", 100)));
        }
    }

    files.CopySummary();
    globalWatch.Stop();

    Console.WriteLine($"Total time: {TimeSpan.FromMilliseconds(globalWatch.ElapsedMilliseconds)}");
}
catch (Exception e)
{
    Console.WriteLine($"Failed with exception {e.Message}. Stack: {e.StackTrace}");
}

return 0;