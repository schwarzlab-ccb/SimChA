using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.Simulation;
using SimChA.IO;
using MathNet.Numerics.Statistics;
using System.Diagnostics;
namespace SimChA.Optimization;
using SimChA.EventData;

public class Optimizer
{   
    private Dictionary<string, List<CopyNumber>> ObservedCNPs { get; }
    private Dictionary<string, List<CopyNumber>> SimulatedCNPs { get; set;}
    protected GenRef GenRef { get; }
    protected readonly Random Rnd;
    protected readonly int Repeats;
    protected readonly SimParams SimParams;
    
    public Optimizer(SimParams simParams, Random rnd, int repeats, GenRef genRef, List<Sample> observedData)
    {
        SimParams = simParams;
        Rnd = rnd;
        Repeats = repeats;
        GenRef = genRef;
        ObservedCNPs = GetCNPs(observedData);
        SimulatedCNPs = new Dictionary<string, List<CopyNumber>>();
    }

    public virtual SimParams Optimize()
        => FindBestParams(5000, 0.01); // 1000 samples, 1% step size

    private double GetScore(Dictionary<string, List<CopyNumber>> cnps)
    {
        var distance = GetEventDistance(cnps);
        return distance;
    }

    private Dictionary<string, List<CopyNumber>> GenerateCNPs(SimParams simParams)
    {
        var samples = GenerateSimulatedData(simParams);
        return GetCNPs(samples);
    }

    private SimParams FindBestParams(int numSamples, double stepFactor)
    {
        var currentParams = SimParams;
        var currentCNPs = GenerateCNPs(currentParams);
        var currentScore = GetScore(currentCNPs);
        for (int i = 0; i < numSamples; i++)
        {
            var proposedParams = GetProposalParams(currentParams, stepFactor);
            var proposedCNPs = GenerateCNPs(proposedParams);
            var proposedScore = GetScore(proposedCNPs);
            var delta = proposedScore - currentScore;
            var prob = Math.Min(1, Math.Exp(-delta));
            if (Rnd.NextDouble() < prob)
            {
                currentParams = proposedParams;
                currentScore = proposedScore;
            }
        }
        return currentParams;
    }

    private SimParams GetProposalParams(SimParams currentParams, double stepFactor)
    {
        if (currentParams.Signatures is null || currentParams.Signatures.Count == 0)
        {
            throw new Exception("Error in Optimizer. No signatures were provided.");
        }
        var events = currentParams.Signatures["CNs"].Events;
        var index = Rnd.Next(events.Count);
        // Modify the relative weight of the event
        var sign = Rnd.NextDouble() < 0.5 ? -1 : 1;
        var newProb = events[index].Prob * (1 + sign * Rnd.NextDouble() * stepFactor);
        while (Math.Abs(newProb - events[index].Prob) <= double.Epsilon)
        {
            sign = Rnd.NextDouble() < 0.5 ? -1 : 1;
            newProb = events[index].Prob * (1 + sign * Rnd.NextDouble() * stepFactor);
        }
        // The other event parameters are multiplicatively modified (with the same relative factor)
        var newFactor = (1 - newProb)/(1 - events[index].Prob);
        var newEvents = events.Select(e => e with { Prob = e.Prob * newFactor }).ToList();
        newEvents[index] = newEvents[index] with { Prob = newProb };
        var newSignature = new Signature(1, newEvents);
        return currentParams with { Signatures = new Dictionary<string, Signature> { ["CNs"] = newSignature } };
    }


    private double GetEventDistance(Dictionary<string, List<CopyNumber>> cnps)
    {
        var segDist = GetSegLengthDistance(cnps);
        var cpDist = GetChangepointDistance(cnps);
        var bpDist = GetBreakpointDistance(cnps);
        var majDist = GetMajMinCNDistance(cnps, true);
        var minDist = GetMajMinCNDistance(cnps, false);
        return (segDist + cpDist + bpDist + majDist + minDist)/5;
    }

    private List<Sample> GenerateSimulatedData(SimParams currentParams)
    {
        if (currentParams.Signatures is null || currentParams.Signatures.Count == 0)
        {
            throw new Exception("No signatures were provided.");
        }
        Validators.ValidateSignatures(currentParams.Signatures);

        var samples = Converters.MakeSamples(Rnd, Repeats, currentParams.EventCount, 
            currentParams.EventDist, currentParams.Signatures, currentParams.Sex, currentParams.MCTarget);
        var simulator = new Simulator(Rnd, GenRef);
        foreach (var sample in samples)
        {
            simulator.SampleEvents(sample);
        }
        return samples;
    }
    private double GetSegLengthDistance(Dictionary<string, List<CopyNumber>> simCNPs)
    {
        var (obsValues, obsMax) = SummaryFeatures.GetSegLengths(ObservedCNPs);
        var (simValues, simMax) = SummaryFeatures.GetSegLengths(simCNPs);
        var histMax = Math.Max(obsMax, simMax);
        var histMin = 0;
        var histBins = 200;
        return CalculateDistance(obsValues, simValues, histBins, histMin, histMax);
    }

    private double GetChangepointDistance(Dictionary<string, List<CopyNumber>> simCNPs)
    {
        var (obsValues, obsMax) = SummaryFeatures.GetChangepointInfo(ObservedCNPs);
        var (simValues, simMax)  = SummaryFeatures.GetChangepointInfo(simCNPs);
        var histMax = Math.Max(obsMax, simMax);
        var histMin = 0;
        var histBins = 50;
        return CalculateDistance(obsValues, simValues, histBins, histMin, histMax);
    }
    private double GetBreakpointDistance(Dictionary<string, List<CopyNumber>> simCNPs)
    {
        var (obsValues, obsMax) = SummaryFeatures.GetBreakpointsPerChromosome(ObservedCNPs);
        var (simValues, simMax)  = SummaryFeatures.GetBreakpointsPerChromosome(simCNPs);
        var histMax = Math.Max(obsMax, simMax);
        var histMin = 0;
        var histBins = 50;
        return CalculateDistance(obsValues, simValues, histBins, histMin, histMax);
    }
    private double GetMajMinCNDistance(Dictionary<string, List<CopyNumber>> simCNPs, bool getMajor)
    {
        var (obsValues, obsMax) = SummaryFeatures.GetMajMinCNs(ObservedCNPs, getMajor);
        var (simValues, simMax) = SummaryFeatures.GetMajMinCNs(simCNPs, getMajor);
        var histMax = Math.Max(obsMax, simMax);
        var histMin = 0;
        var histBins = 50;
        return CalculateDistance(obsValues, simValues, histBins, histMin, histMax);
    }
    protected static double CalculateDistance(List<double> data, List<double> sim, int bins, int min, double max)
    {
        var dataHist = new Histogram(data, bins, min, max);
        var simHist  = new Histogram(sim, bins, min, max);
        return StatisticMeasures.WassersteinDistance(dataHist, simHist);
    }
    protected Dictionary<string, List<CopyNumber>> GetCNPs(List<Sample> samples)
    {
        var cnps = new Dictionary<string, List<CopyNumber>>();
        foreach (var sample in samples)
        {
            foreach (var clone in sample.Clones)
            {
                var cnp = CopyNumbers.CalcCopyNumbers(GenRef, sample.Kars[clone.CloneId], sample.Kars[clone.CloneId].SexXX).ToList();
                cnps[sample.SampleId] = cnp;
            }
        }
        return cnps;
    }
}