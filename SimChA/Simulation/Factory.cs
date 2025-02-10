using SimChA.Computation;
using SimChA.Data;
using SimChA.EventData;
using SimChA.IO;

namespace SimChA.Simulation;

public static class Factory
{
    public static Simulator GetSimulator(Random rnd, GenRef genRef, SimChAConfig simChAConfig, SelectionMode selMode)
    {
        var sampleParams = simChAConfig.SimParams;
        var fitParams = simChAConfig.FitParams;
        switch (selMode)
        {
            case SelectionMode.MetropolisHastings:
                var mcParams = simChAConfig.MHParams ?? throw new Exception(
                    "Error: MCParams not set. Cannot perform MC sampling. Please set MCParams in the config file.");
                return new MHSimulator(rnd, genRef, sampleParams, fitParams, mcParams);

            case SelectionMode.SimulatedAnnealing:
                var evoParams = simChAConfig.EvoParams ?? throw new Exception(
                    "Error: EvoParams not set. Cannot perform evolution without evolution parameters. Please set EvoParams in the config file.");
                return new SASimulator(rnd, genRef, sampleParams, fitParams, evoParams);

            case SelectionMode.MonteCarlo:
                return new Simulator(rnd, genRef, sampleParams, fitParams);

            default:
                throw new Exception("Error: No selection mode set.");
        }
    }

    public static List<Sample> ReadSamples(Random rnd, GenRef genRef, SimChAConfig config, CmdOptions options)
    {
        var sampleParams = config.SimParams;
        if (config.SimParams == null)
        {
            throw new Exception("No simulation parameters found. Please set \"SimParams\" in the config JSON.");
        }

        if (config.FitParams == null)
        {
            throw new Exception("No fitness parameters found. Please set \"FitParams\" in the config JSON.");
        }
        
        if (options.ExecMode == ExecMode.Profiles)
        {
            Console.WriteLine("Reading samples from data:");
            var profiles = FileIO.ReadProfiles(genRef, options.CNProfiles, sampleParams.AutosomesOnly);
            return Converters.SamplesFromProfiles(profiles);
        }
        
        var validSigs = ValidateSignatures(config.Signatures);
        if (options.ExecMode == ExecMode.Repeats)
        {
            Console.WriteLine("Creating random samples:");
            return Converters.MakeSamples(rnd, options.Repeats, validSigs, sampleParams.Sex, sampleParams.AutosomesOnly);
        }
        if (options.ExecMode == ExecMode.Tree)
        {
            Console.WriteLine("Reading samples from a clone file:");
            var inClones = FileIO.ReadCloneTree(options.CloneTreeFile, options.MHMode);
            var (cnEventPs, mixture) = Converters.PropagateSigs(validSigs);
            string sampleName = Path.GetFileNameWithoutExtension(options.CloneTreeFile);
            var sex = sampleParams.AutosomesOnly ? SexType.Any : Sampling.GetSex(rnd, sampleParams.Sex);
            var treeSample = new Sample(sampleName, sex, inClones, cnEventPs, mixture);
            return new List<Sample> { treeSample };
        }
        throw new Exception("Unknown execution mode.");
    }

    private static void ValidateEvent(CNEventPars cnEventPars)
    {
        switch (cnEventPars.Type)
        {
            case CNEventType.ChromDeletion:
            case CNEventType.ChromDuplication:
            case CNEventType.Pass:
            case CNEventType.WholeGenomeDoubling:
            case CNEventType.BreakageFusionBridge:
            case CNEventType.SNV:
            case CNEventType.ArmDeletion:
            case CNEventType.ArmDuplication:
                break;
                
            case CNEventType.InternalDeletion:
            case CNEventType.InternalDuplication:
            case CNEventType.InternalInversion:
            case CNEventType.InvertedDuplication:
            case CNEventType.Translocation:
            case CNEventType.Chromothripsis:
            case CNEventType.Chromoplexy:
            case CNEventType.CentromereBoundDeletion:
            case CNEventType.CentromereBoundDuplication:
            case CNEventType.TailDeletion:
            case CNEventType.TailDuplication:
                if (cnEventPars.Size <= 0) 
                    throw new Exception($"Event {cnEventPars.Type} does not have a Size parameter. E.g. \"Size\": 1000000");
                break;
            
            case CNEventType.TIChain:
            case CNEventType.TIBridge:
            case CNEventType.TICycle:
            case CNEventType.Pyrgo:
            case CNEventType.Rigma:
                if (cnEventPars.Size <= 0) 
                    throw new Exception($"Event {cnEventPars.Type} does not have a Size parameter. E.g. \"Size\": 1000000");
                if (cnEventPars.Frag <= 0) 
                    throw new Exception($"Event {cnEventPars.Type} does not have a Size parameter. E.g. \"Frag\": 5");
                break;
            
            default:
                throw new ArgumentOutOfRangeException($"Unknown event type {cnEventPars.Type}");
        }
    }

    private static List<Signature> ValidateSignatures(List<Signature>? signatures)
    {
        if (signatures is null)
        {
            throw new Exception("No signatures found.");
        }
        var validSigs = new List<Signature>();
        
        foreach (var sig in signatures.Where(sig => sig.Prob > 0 && sig.Events.Any(e => e.Prob > 0)))
        {
            if (sig.Events is null || sig.Events.Count == 0)
            {
                throw new Exception($"Signature {sig.Name} does not have any events.");
            }
            double probSum = sig.Events.Sum(e => e.Prob);
            if (probSum <= 0)
            {
                throw new Exception($"Signature {sig.Name} has a total probability of {probSum}.");
            }
            foreach(var cnEventP in sig.Events)
            {
                ValidateEvent(cnEventP);
            }
            validSigs.Add(sig);
        }
        if (signatures.Count == 0)
        {
            throw new Exception("No valid signatures found.");
        }
        return validSigs;
    }
}
