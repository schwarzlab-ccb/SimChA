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
    protected GenRef GenRef { get; }
    protected readonly Random Rnd;
    protected readonly int Repeats;
    protected readonly SimParams SimParams;
    protected readonly OptimizationParams OptimizationParams;
    
    public Optimizer(SimParams simParams, Random rnd, int repeats, GenRef genRef, List<Sample> observedData)
    {
        SimParams = simParams;
        Rnd = rnd;
        Repeats = repeats;
        GenRef = genRef;
        ObservedCNPs = GetCNPs(observedData);
        OptimizationParams = SimParams.OptimizationParams ?? throw new Exception("Error in Optimizer. OptimizationParams not set.");
    }

    public virtual SimParams Optimize(FileIO files)
        => FindBestParams(files);

    private double GetScore(Dictionary<string, List<CopyNumber>> cnps)
    {
        var segDist = GetSegLengthDistance(cnps);
        var cpDist = GetChangepointDistance(cnps);
        var bpDist = GetBreakpointDistance(cnps);
        var majDist = GetMajMinCNDistance(cnps, true);
        var minDist = GetMajMinCNDistance(cnps, false);
        return segDist;//(segDist + cpDist + bpDist + majDist + minDist)/5;
    }

    private Dictionary<string, List<CopyNumber>> GenerateCNPs(SimParams currentParams)
    {
        var samples = GenerateSimulatedData(currentParams);
        return GetCNPs(samples);
    }

    private SimParams FindBestParams(FileIO files)
    {
        var currentParams = GetProposalParams(SimParams);
        var currentCNPs = GenerateCNPs(currentParams);
        var currentScore = GetScore(currentCNPs);
        int counter = 0;
        for (int i = 0; i < OptimizationParams.NumSamplesTotal; i++)
        {
            Console.WriteLine($"Iteration {i+1} of {OptimizationParams.NumSamplesTotal}");
            var proposedParams = GetProposalParams(currentParams);
            var proposedCNPs = GenerateCNPs(proposedParams);
            var proposedScore = GetScore(proposedCNPs);
            var delta = OptimizationParams.AcceptanceFactor*(proposedScore - currentScore)/currentScore;
            var prob = Math.Min(1, Math.Exp(-delta));
            if (prob < 1)
            {
                Console.WriteLine("hello");
            }
            if (Rnd.NextDouble() < prob)
            {
                currentParams = proposedParams;
                currentScore = proposedScore;
                counter++;
                if (OptimizationParams.WriteIntermediate && counter % OptimizationParams.WriteFrequency == 0)
                    files.WriteSimParams(currentParams, $"params_{counter}.json");
            }
        }
        return currentParams;
    }

    private SimParams GetProposalParams(SimParams currentParams)
    {
        if (currentParams.Signatures is null || currentParams.Signatures.Count == 0)
        {
            throw new Exception("Error in Optimizer. No signatures were provided.");
        }
        var events = currentParams.Signatures["CNs"].Events.Where(e => e.Prob > 0).ToList();
        var totalWeight = events.Sum(e => e.Prob);
        // Choose an event to modify
        // We add an extra two possible indices to account for the InternalDuplication and
        // InternalDeletion Length Parameters
        var index = Rnd.Next(events.Count + 2);
        var sign = Rnd.NextDouble() < 0.5 ? -1 : 1;
        var nTries = 0;
        // Modify internal length scales
        if (index >= events.Count)
        {
            long oldSize = 0;
            if (index == events.Count)
            {
                var internalDel = events.Find(e => e.Type == CNEventType.InternalDeletion) ?? throw new Exception("Error in Optimizer. No internal deletion event found.");
                // Get the actual index
                index = events.IndexOf(internalDel);
                oldSize = internalDel.Size;
            }
            else if (index == events.Count + 1)
            {
                var internalDup = events.Find(e => e.Type == CNEventType.InternalDuplication) ?? throw new Exception("Error in Optimizer. No internal duplication event found.");
                index = events.IndexOf(internalDup);
                oldSize = internalDup.Size;
            }
            else
            {
                throw new Exception("Error in Optimizer. Index out of range.");
            }
            // Modify the internal duplication/deletion length parameter
            var newSize = (long) (oldSize * (1 + sign * Rnd.NextDouble() * OptimizationParams.StepFactor));
            while (newSize < 0 && Math.Abs(newSize-oldSize)/oldSize <= double.Epsilon && nTries < OptimizationParams.MaxTries)
            {
                nTries++;
                sign = Rnd.NextDouble() < 0.5 ? -1 : 1;
                newSize = (long) (oldSize * (1 + sign * Rnd.NextDouble() * OptimizationParams.StepFactor));
            }
            var newEvents = new List<CNEventPars>(events);
            newEvents[index] = newEvents[index] with { Size = newSize };
            var newSignature = new Signature(1, newEvents);
            return currentParams with { Signatures = new Dictionary<string, Signature> { ["CNs"] = newSignature } };
        }
        else
        {
            // Modify the relative weight of the event
            var newProb = events[index].Prob * (1 + sign * Rnd.NextDouble() * OptimizationParams.StepFactor);

            while (Math.Abs(newProb - events[index].Prob)/events[index].Prob <= double.Epsilon && nTries < OptimizationParams.MaxTries)
            {
                nTries++;
                sign = Rnd.NextDouble() < 0.5 ? -1 : 1;
                newProb = events[index].Prob * (1 + sign * Rnd.NextDouble() * OptimizationParams.StepFactor);
            }
            // The other event parameters are multiplicatively modified (with the same relative factor)
            var newFactor = (totalWeight - newProb)/(totalWeight - events[index].Prob);
            var newEvents = events.Select(e => e with { Prob = e.Prob * newFactor }).ToList();
            newEvents[index] = newEvents[index] with { Prob = newProb };
            var newSignature = new Signature(1, newEvents);
            return currentParams with { Signatures = new Dictionary<string, Signature> { ["CNs"] = newSignature } };
        }
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
        var histBins = histMax;
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