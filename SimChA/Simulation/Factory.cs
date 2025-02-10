using SimChA.Computation;
using SimChA.Data;
using SimChA.EventData;
using SimChA.IO;

namespace SimChA.Simulation;

public static class Factory
{
    public static Simulator GetSimulator(Random rnd, GenRef genRef, SimParams simParams, SelectionMode selMode)
    {
        var sampleParams = simParams.SampleParams;
        var fitParams = simParams.FitParams;
        switch (selMode)
        {
            case SelectionMode.MetropolisHastings:
                var mcParams = simParams.MHParams ?? throw new Exception(
                    "Error: MCParams not set. Cannot perform MC sampling. Please set MCParams in the config file.");
                return new MHSimulator(rnd, genRef, sampleParams, fitParams, mcParams);

            case SelectionMode.SimulatedAnnealing:
                var evoParams = simParams.EvoParams ?? throw new Exception(
                    "Error: EvoParams not set. Cannot perform evolution without evolution parameters. Please set EvoParams in the config file.");
                return new SASimulator(rnd, genRef, sampleParams, fitParams, evoParams);

            case SelectionMode.MonteCarlo:
                return new Simulator(rnd, genRef, sampleParams, fitParams);

            default:
                throw new Exception("Error: No selection mode set.");
        }
    }

    public static List<Sample> ReadSamples(Random rnd, GenRef genRef, SimParams simParams, CmdOptions options)
    {
        var sampleParams = simParams.SampleParams;
        if (options.ExecMode == ExecMode.Profiles)
        {
            Console.WriteLine("Reading samples from data:");
            var profiles = FileIO.ReadProfiles(genRef, options.CNProfiles, sampleParams.AutosomesOnly);
            return Converters.SamplesFromProfiles(profiles);
        }
        
        var validSigs = ValidateSignatures(simParams.Signatures);

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
            var sex = sampleParams.AutosomesOnly ? SexType.None : Sampling.GetSex(rnd, sampleParams.Sex);
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

    private static Dictionary<string, Signature> ValidateSignatures(Dictionary<string, Signature>? signatures)
    {
        if (signatures is null)
        {
            throw new Exception("No signatures found.");
        }
        var validSigs = new Dictionary<string, Signature>();
        
        foreach (var sig in signatures.Where(sig => sig.Value.Prob > 0 && sig.Value.Events.Any(e => e.Prob > 0)))
        {
            if (sig.Value.Events is null || sig.Value.Events.Count == 0)
            {
                throw new Exception($"Signature {sig.Key} does not have any events.");
            }
            double probSum = sig.Value.Events.Sum(e => e.Prob);
            if (probSum <= 0)
            {
                throw new Exception($"Signature {sig.Key} has a total probability of {probSum}.");
            }
            foreach(var cnEventP in sig.Value.Events)
            {
                ValidateEvent(cnEventP);
            }
            validSigs[sig.Key] = sig.Value;
        }
        if (signatures.Count == 0)
        {
            throw new Exception("No valid signatures found.");
        }
        return validSigs;
    }
}
