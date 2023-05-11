// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.EventData;
using SimChA.Misc;
using SimChA.Simulation;

namespace SimChA.IO;

public static class Converters
{
    // This will make events conditional on the probability of the signature
    // If probs are not null, they are used instead of the signature probs and have to have the same length
    public static List<CNEventPars> PropagateSigs(List<Signature> signatures, List<double>? probs = null)
    {
        var events = new List<CNEventPars>();
        for (var i = 0; i < signatures.Count; i++)
        {
            if (probs != null)
            {
                if (probs.Count != signatures.Count)
                {
                    throw new ArgumentException("probs must have the same length as signatures");
                }
                if (Math.Abs(probs.Sum() - 1) > double.Epsilon * probs.Count)
                {
                    throw new ArgumentException("probs must sum to 1");
                }
            }
            var signature = signatures[i];
            double sigProbSum = signatures.Sum(sig => sig.Prob);
            double sigProb = probs == null ? signature.Prob / sigProbSum : probs[i];
            
            var selectedEvs = signature.Events.Where(ev => ev.Prob > 0).ToList();
            double evsProbSum = selectedEvs.Sum(ev => ev.Prob);
            events.AddRange(selectedEvs.Select(cnEventP => cnEventP with
            {
                Prob = cnEventP.Prob / evsProbSum * sigProb
            }));
        }
        return events;    
    }
    
    public static List<Sample> MakeSamples(Random rnd, int repeats, int meanDist, Distribution distribution, List<Signature> sigs)
    {
        var samples = new List<Sample>();
        var selectedSigs = sigs.Where(s => s.Prob > 0).ToList();
        double[] sigProbs = sigs.Select(s => s.Prob).ToArray();
        var mixture = Sampling.CreateRandomMixture(rnd, sigProbs);
        for (var i = 0; i < repeats; i++)
        {
            double dist = Sampling.SampleDist(rnd, distribution);
            var mutCount = (int)Math.Round(meanDist * dist);
            var clone = new CloneIn(0, -1, mutCount, 1); // TODO: Specify fitness target
            var events = PropagateSigs(selectedSigs, mixture);
            bool sexXX = rnd.CoinFlip();
            var sample = new Sample($"Sample{i + 1}", sexXX, new List<CloneIn> { clone }, events);
            samples.Add(sample);
        }
        return samples;
    }
}