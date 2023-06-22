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
    public static (List<CNEventPars>, Dictionary<string, double> mixture) PropagateSigs(List<Signature> signatures, List<double>? probs = null)
    {
        var events = new List<CNEventPars>();
        var mixture = new Dictionary<string, double>();
        for (int i = 0; i < signatures.Count; i++)
        {
            if (probs != null)
            {
                if (probs.Count != signatures.Count)
                {
                    throw new ArgumentException("probs must have the same length as signatures");
                }
                if (Math.Abs(probs.Sum() - 1) > 1e-9)
                {
                    throw new ArgumentException("probs must sum to 1");
                }
            }
            var signature = signatures[i];
            double sigProbSum = signatures.Sum(sig => sig.Prob);
            double sigProb = probs == null ? signature.Prob / sigProbSum : probs[i];
            mixture.Add(signature.Id, sigProb);
            
            var selectedEvs = signature.Events.Where(ev => ev.Prob > 0).ToList();
            double evsProbSum = selectedEvs.Sum(ev => ev.Prob);
            events.AddRange(selectedEvs.Select(cnEventP => cnEventP with
            {
                Prob = cnEventP.Prob / evsProbSum * sigProb
            }));
        }
        return (events, mixture);    
    }
    
    public static List<Sample> MakeSamples(Random rnd, int repeats, int meanDist, Distribution distribution, List<Signature> sigs, SexEnum sex)
    {
        var samples = new List<Sample>();
        var selectedSigs = sigs.Where(s => s.Prob > 0).ToList();
        string[] sigNames = sigs.Select(s => s.Id).ToArray();
        double[] sigProbs = sigs.Select(s => s.Prob).ToArray();
        for (int i = 0; i < repeats; i++)
        {
            double dist = Sampling.SampleDist(rnd, distribution);
            int mutCount = (int) Math.Round(meanDist * dist);
            double fitnessTarget = Sampling.SampleDist(rnd, Distribution.Exponential) * 0.8 + 1;
            var clone = new CloneIn(0, -1, mutCount, fitnessTarget); 
            var dirichlet = Sampling.CreateRandomMixture(rnd, sigProbs);
            var (events, mixture) = PropagateSigs(selectedSigs, dirichlet);
            var sample = new Sample($"sample_{i + 1}", Sampling.GetBinarySex(rnd, sex), new List<CloneIn> { clone }, events, mixture);
            samples.Add(sample);
        }
        return samples;
    }
}