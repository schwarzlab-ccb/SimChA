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

    public static List<Sample> MakeSamples(
        Random rnd,
        int repeats,
        List<Signature> sigs,
        SexType sex,
        bool autosomesOnly)
    {
        var samples = new List<Sample>();
        string[] sigNames = sigs.Select(s => s.Name).ToArray();
        double[] sigProbs = sigs.Select(s => s.Prob).ToArray();
        for (int i = 0; i < repeats; i++)
        {
            // Sample to have at least 1 event
            var clone = new CloneIn($"clone_{i + 1}", $"clone_{i + 1}", -1, -1);
            var dirichlet = Sampling.CreateRandomMixture(rnd, sigProbs);
            var namedProbs = sigNames.Zip(dirichlet).ToDictionary(s => s.First, s => s.Second);
            var (events, mixture) = PropagateSigs(sigs, namedProbs);
            var sampleSex = autosomesOnly ? SexType.Any : Sampling.GetSex(rnd, sex);
            var sample = new Sample($"sample_{i + 1}", sampleSex, new List<CloneIn> { clone }, events, mixture);
            samples.Add(sample);
        }
        return samples;
    }

    public static List<Sample> SamplesFromProfiles(Dictionary<string, Karyotype> profiles)
        => (from profile in profiles
            let clones = new List<CloneIn> { new("0", "-1", 0, 0) }
            select new Sample(profile.Key, profile.Value.Sex, clones)
                { Kars = { ["0"] = profile.Value } }).ToList();

}