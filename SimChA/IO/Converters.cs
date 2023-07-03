// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;
using SimChA.Simulation;

namespace SimChA.IO;

public static class Converters
{
    // This will make events conditional on the probability of the signature
    // If probs are not null, they are used instead of the signature probs and have to have the same length
    public static (List<CNEventPars>, Dictionary<string, double> mixture) PropagateSigs(Dictionary<string, Signature> signatures, Dictionary<string, double>? probs = null)
    {
        var events = new List<CNEventPars>();
        var mixture = new Dictionary<string, double>();
        double sigProbSum = signatures.Sum(sig => sig.Value.Prob);
        foreach(string sigId in signatures.Keys) 
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
            var signature = signatures[sigId];
            double sigProb = probs == null ? signature.Prob / sigProbSum : probs[sigId];
            mixture.Add(sigId, sigProb);
            
            var selectedEvs = signature.Events.Where(ev => ev.Prob > 0).ToList();
            double evsProbSum = selectedEvs.Sum(ev => ev.Prob);
            events.AddRange(selectedEvs.Select(cnEventP => cnEventP with
            {
                Prob = cnEventP.Prob / evsProbSum * sigProb
            }));
        }
        return (events, mixture);    
    }
    
    public static List<Sample> MakeSamples(Random rnd, int repeats, int meanDist, Distribution distribution, Dictionary<string, Signature> sigs, SexEnum sex)
    {
        var samples = new List<Sample>();
        var selectedSigs = sigs.Where(s => s.Value.Prob > 0).ToDictionary(s => s.Key, s => s.Value);
        string[] sigNames = sigs.Select(s => s.Key).ToArray();
        double[] sigProbs = sigs.Select(s => s.Value.Prob).ToArray();
        for (int i = 0; i < repeats; i++)
        {
            double dist = Sampling.SampleDist(rnd, distribution);
            int mutCount = (int) Math.Round(meanDist * dist);
            double fitnessTarget = Sampling.SampleDist(rnd, Distribution.Exponential) * 0.8 + 1;
            var clone = new CloneIn(0, -1, mutCount, fitnessTarget); 
            var dirichlet = Sampling.CreateRandomMixture(rnd, sigProbs);
            var namedProbs = sigNames.Zip(dirichlet).ToDictionary(s => s.First, s => s.Second);
            var (events, mixture) = PropagateSigs(selectedSigs, namedProbs);
            var sample = new Sample($"sample_{i + 1}", Sampling.GetBinarySex(rnd, sex), new List<CloneIn> { clone }, events, mixture);
            samples.Add(sample);
        }
        return samples;
    }
}