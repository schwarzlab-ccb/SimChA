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
                if (cnEventPars.Frac <= 0) 
                    throw new Exception($"Event {cnEventPars.Type} does not have a Size parameter. E.g. \"Size\": 1000000");
                break;
            
            case CNEventType.TIChain:
            case CNEventType.TIBridge:
            case CNEventType.TICycle:
            case CNEventType.Pyrgo:
            case CNEventType.Rigma:
                if (cnEventPars.Frac <= 0) 
                    throw new Exception($"Event {cnEventPars.Type} does not have a Size parameter. E.g. \"Size\": 1000000");
                if (cnEventPars.Frag <= 0) 
                    throw new Exception($"Event {cnEventPars.Type} does not have a Size parameter. E.g. \"Frag\": 5");
                break;
            
            default:
                throw new ArgumentOutOfRangeException($"Unknown event type {cnEventPars.Type}");
        }
    }

    public static List<Signature> ValidateSignatures(List<Signature>? signatures)
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

    public static (List<CNEventPars>, Dictionary<string, double> mixture) MixSignatures(List<Signature> signatures, Dictionary<string, double>? probs = null)
    {
        if (probs != null)
        {
            if (probs.Count != signatures.Count)
            {
                throw new ArgumentException("probs must have the same length as signatures");
            }
            if (Math.Abs(probs.Values.Sum() - 1) > 1e-9)
            {
                throw new ArgumentException("probs must sum to 1");
            }
        }
        var events = new List<CNEventPars>();
        var mixture = new Dictionary<string, double>();
        double sigProbSum = signatures.Sum(sig => sig.Prob);
        foreach((string sigId, double prob, var cnEventParsList) in signatures)
        {
            double sigProb = probs == null ? prob / sigProbSum : probs[sigId];
            mixture.Add(sigId, sigProb);
            
            var selectedEvs = cnEventParsList.Where(ev => ev.Prob > 0).ToList();
            double evsProbSum = selectedEvs.Sum(ev => ev.Prob);
            events.AddRange(selectedEvs.Select(cnEventP => cnEventP with
            {
                Prob = cnEventP.Prob / evsProbSum * sigProb
            }));
        }
        return (events, mixture);    
    }

    public static List<CNEventPars> NormalizeEvents(List<CNEventPars> events)
    {
        double probSum = events.Sum(ev => ev.Prob);
        return events.Select(ev => ev with { Prob = ev.Prob/probSum }).ToList();
    }
}
