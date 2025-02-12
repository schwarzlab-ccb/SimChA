using SimChA.Data;
using SimChA.EventData;
using SimChA.Simulation;

namespace SimChA.Computation;

public static class Converters
{
    // This will make events conditional on the probability of the signature
    // If probs are not null, they are used instead of the signature probs and have to have the same length
    public static (List<CNEventPars>, Dictionary<string, double> mixture) PropagateSigs(List<Signature> signatures, Dictionary<string, double>? probs = null)
    {
        var events = new List<CNEventPars>();
        var mixture = new Dictionary<string, double>();
        double sigProbSum = signatures.Sum(sig => sig.Prob);
        foreach(var sig in signatures)
        {
            string sigId = sig.Name;
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
            double sigProb = probs == null ? sig.Prob / sigProbSum : probs[sigId];
            mixture.Add(sigId, sigProb);
            
            var selectedEvs = sig.Events.Where(ev => ev.Prob > 0).ToList();
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