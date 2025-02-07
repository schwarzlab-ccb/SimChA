using SimChA.DataTypes;
using SimChA.IO;

namespace SimChA.Simulation;

public static class Factory
{
    public static Simulator GetSimulator(Random rnd, GenRef genRef, SimParams simParams, SelectionMode selMode)
    {
        var fitParams = simParams.Fitness;
        switch (selMode)
        {
            case SelectionMode.MetroplisHastings:
                var mcParams = simParams.MCParams ?? throw new Exception(
                    "Error: MCParams not set. Cannot perform MC sampling. Please set MCParams in the config file.");
                return new MHSimulator(rnd, genRef, fitParams, mcParams);

            case SelectionMode.SimulatedAnnealing:
                var evoParams = simParams.EvoParams ?? throw new Exception(
                    "Error: EvoParams not set. Cannot perform evolution without evolution parameters. Please set EvoParams in the config file.");
                return new SASimulator(rnd, genRef, fitParams, evoParams);

            case SelectionMode.MonteCarlo:
                return new Simulator(rnd, genRef, fitParams);

            default:
                throw new Exception("Error: No selection mode set.");
        }
    }

    public static List<Sample> ReadSamples(Random rnd, GenRef genRef, SimParams simParams, CmdOptions options)
    {
        switch (options.ExecMode)
        {
            case ExecMode.Profiles:
                Console.WriteLine("Reading samples from data:");
                var profiles = FileIO.ReadProfiles(genRef, options.CNProfiles, simParams.AutosomesOnly);
                // TODO: Needs to implement autosomes only
                return Converters.SamplesFromProfiles(profiles);

            case ExecMode.Repeats:
                Console.WriteLine("Creating random samples:");
                var sigs = Validators.ValidateSignatures(simParams.Signatures);
                return Converters.MakeSamples(rnd, options.Repeats, sigs, simParams.Sex, simParams.AutosomesOnly);

            case ExecMode.Tree:
                Console.WriteLine("Reading samples from a clone file:");
                var treeSigs = Validators.ValidateSignatures(simParams.Signatures);
                var inClones = FileIO.ReadCloneTree(options.CloneTreeFile, options.MHMode);
                var (cnEventPs, mixture) = Converters.PropagateSigs(treeSigs);
                string sampleName = Path.GetFileNameWithoutExtension(options.CloneTreeFile);
                var sex = simParams.AutosomesOnly ? SexEnum.None : Sampling.GetSex(rnd, simParams.Sex);
                var treeSample = new Sample(sampleName, sex, inClones, cnEventPs, mixture);
                return new List<Sample> { treeSample };

            default:
                throw new Exception("Error: No execution mode set.");
        }
    }
}
