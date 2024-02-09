using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.Simulation;
using SimChA.IO;
using MathNet.Numerics.Statistics;
namespace SimChA.Optimization;
using SimChA.EventData;
using System;

public class Optimizer
{   
    protected Dictionary<string, List<CopyNumber>> ObservedCNPs { get; }
    protected GenRef GenRef { get; }
    protected readonly Random Rnd;
    protected readonly int Repeats;
    protected readonly SimParams SimParams;
    protected readonly OptimizationParams OptimizationParams;
    protected Dictionary<string, bool> IsFemaleObservedDict {get; set;}
    protected Dictionary<string, bool> IsFemaleSimulatedDict;
    private bool IncludeSexChromosomes { get; }
    private bool BreakpointsPerChrom {get;}
    private long BPBinSize {get;}
    
    public Optimizer(SimParams simParams, Random rnd, int repeats, GenRef genRef, List<Sample> observedData, bool includeSexChromosomes)
    {
        SimParams = simParams;
        Rnd = rnd;
        Repeats = repeats;
        GenRef = genRef;
        ObservedCNPs = GetCNPs(observedData);
        IsFemaleObservedDict = observedData.ToDictionary(s => s.SampleId, s => s.SexXX);
        OptimizationParams = SimParams.OptimizationParams ?? throw new Exception("Error in Optimizer. OptimizationParams not set.");
        IncludeSexChromosomes = includeSexChromosomes;
        BreakpointsPerChrom = OptimizationParams.BreakpointsPerChrom;
        BPBinSize = OptimizationParams.BreakpointsBinSize;
        if (SimParams.Signatures is null || SimParams.Signatures.Count == 0)
        {
            throw new Exception("Error in Optimizer. No signatures were provided.");
        }
    }

    public virtual SimParams Optimize(FileIO files)
        => FindBestParams(files);
    
    public double GetABCDistance()
    {
        var simCNPs = GenerateCNPs(SimParams);
        return GetScore(simCNPs);
    }

    public double GetScore(Dictionary<string, List<CopyNumber>> cnps)
    {
        var totalDist = 0.0;
        if (OptimizationParams.UseMeanSeg)
        {
            var segDist = GetMeanSegDistance(cnps);
            totalDist += segDist*segDist;
        }
        if (OptimizationParams.UsePloidy)
        {
            var ploidyDist = GetPloidyDistance(cnps);
            totalDist += ploidyDist*ploidyDist;
        }
        if (OptimizationParams.UseBreakpoints)
        {
            var bpDist = GetBreakpointDistance(cnps);
            totalDist += bpDist*bpDist;
        }
        //var copyNumberMatrix = SummaryFeatures.GetChrCopyNumberMatrix(GenRef.AllChrs, cnps);
        //var mkv = SummaryFeatures.GetMKV(copyNumberMatrix);
        //var aneuploidy = SummaryFeatures.GetAverageAneuploidy(copyNumberMatrix);
        return totalDist;
    }

    private Dictionary<string, List<CopyNumber>> GenerateCNPs(SimParams currentParams)
    {
        var samples = GenerateSimulatedData(currentParams);
        IsFemaleSimulatedDict = samples.ToDictionary(s => s.SampleId, s => s.SexXX);
        return GetCNPs(samples);
    }

    private SimParams FindBestParams(FileIO files)
    {
        var currentParams = SimParams;//GetProposalParams(SimParams);
        if (OptimizationParams.ResetSeed)
        {
            currentParams = currentParams with {Seed = -1};
        }
        var currentCNPs = GenerateCNPs(currentParams);
        var currentScore = GetScore(currentCNPs);
        var bestParams = currentParams;
        var bestScore = currentScore;
        var counter = 0;
        for (int i = 0; i < OptimizationParams.NumSamplesTotal; i++)
        {
            Console.WriteLine($"Iteration {i+1} of {OptimizationParams.NumSamplesTotal}");
            var proposedParams = GetProposalParams(currentParams);
            var proposedCNPs = GenerateCNPs(proposedParams);
            var proposedScore = GetScore(proposedCNPs);
            var delta = OptimizationParams.AcceptanceFactor*(proposedScore - currentScore)/currentScore;
            var prob = Math.Min(1, Math.Exp(-delta));
            if (Rnd.NextDouble() < prob)
            {
                currentParams = proposedParams;
                currentScore = proposedScore;
            }
            if (proposedScore < bestScore)
            {
                bestParams = proposedParams;
                bestScore = proposedScore;
                files.WriteSimParams(bestParams, $"best_params_{counter}.json");
                counter++;
            }
            if (OptimizationParams.WriteIntermediate && i % OptimizationParams.WriteFrequency == 0)
            {
                files.WriteSimParams(currentParams, $"params_{i}.json");
            }
            if (i % OptimizationParams.GCInterval == 0)
            {
                GC.Collect();
            }
        }
        return bestParams;
    }

    private SimParams GetProposalParams(SimParams currentParams)
    {
        return OptimizationParams.ParamVariationMode switch
        {
            0 => GetAllNewParams(currentParams),
            _ => GetNNewParams(currentParams),
        };
    }

    private SimParams GetNNewParams(SimParams currentParams)
    {
        var n = OptimizationParams.ParamVariationMode;
        var events = currentParams.Signatures["CNs"].Events.Where(e => e.Prob > 0).ToList();
        var newProbs = events.Select(e => e.Prob).ToList();
        var targetWeight = newProbs.Sum();
        var nParams = OptimizationParams.EventWeightsOnly ? events.Count : events.Count + 2;
        if (n > nParams)
        {
            throw new Exception("Error in Optimizer. n is greater than the number of possible parameters.");
        }
        // Select the two event properties to modify
        var indices = Enumerable.Range(0, nParams).OrderBy(x => Rnd.Next()).Take(n).ToList();
        var newEvents = new List<CNEventPars>(events);
        bool reweightNeeded = false;
        var sumNewWeights = 0.0;
        foreach (var index in indices)
        {
            if (index >= events.Count)
            {
                var eventIndex = index == events.Count ? events.IndexOf(events.Find(e => e.Type == CNEventType.InternalDeletion)) : events.IndexOf(events.Find(e => e.Type == CNEventType.InternalDuplication));
                newEvents[eventIndex] = GetNewEventLength(newEvents[eventIndex]);
            }
            else
            {   
                newProbs[index] = GetNewWeight(events[index].Prob);
                sumNewWeights += newProbs[index];
                reweightNeeded = true;
            }
        }
        if (reweightNeeded)
        {
            // Reweight the ones that weren't changed by a common factor
            var newTotal = newProbs.Sum();
            var factor = (targetWeight - sumNewWeights)/(newTotal - sumNewWeights);
            // If the factor is negative, reweight all events
            if (factor <= 0.0)
            {
                newProbs = newProbs.Select(x => targetWeight * x / newTotal).ToList();
            }
            else
            {
                newProbs = newProbs.Select((x,i) => indices.Contains(i) ? x : factor * x).ToList();
                if (newProbs.Any(x => x <= 0.0))
                {
                    throw new Exception("Error in Optimizer. Negative probability.");
                }
            }
        }
        newEvents = newEvents.Select((e, i) => e with { Prob = newProbs[i] }).ToList();
        var newSignature = new Signature(1, newEvents);
        return currentParams with { Signatures = new Dictionary<string, Signature> { ["CNs"] = newSignature } };
    }

    private double GetNewWeight(double oldProb)
    {
        var nTries = 0;
        var newProb = oldProb * (1 + (Rnd.NextDouble() < 0.5 ? -1 : 1) * Rnd.NextDouble() * OptimizationParams.StepFactor);
        while (newProb <= 0.0 && nTries < OptimizationParams.MaxTries)
        {
            nTries++;
            newProb = oldProb * (1 + (Rnd.NextDouble() < 0.5 ? -1 : 1) * Rnd.NextDouble() * OptimizationParams.StepFactor);
        }
        return newProb;
    }

    private long GetNewLength(long oldValue)
        => (long)GetNewWeight(oldValue);
    
    private SimParams GetAllNewParams(SimParams currentParams)
    {
        var events = currentParams.Signatures["CNs"].Events.Where(e => e.Prob > 0).ToList();
        var targetWeight = events.Sum(e => e.Prob);
        var newProbs = new List<double>();
        foreach (var ev in events)
        {
            var weight = GetNewWeight(ev.Prob);
            newProbs.Add(weight);
        }
        var currentTotal = newProbs.Sum();
        newProbs  = newProbs.Select(x => targetWeight * x / currentTotal).ToList();
        var newEvents = events.Select((e, i) => e with { Prob = newProbs[i] }).ToList();

        if (!OptimizationParams.EventWeightsOnly)
        {
            var internalDel = events.Find(e => e.Type == CNEventType.InternalDeletion) ?? throw new Exception("Error in Optimizer. No internal deletion event found.");
            var index = events.IndexOf(internalDel);
            newEvents[index] = GetNewEventLength(newEvents[index]);

            var internalDup = events.Find(e => e.Type == CNEventType.InternalDuplication) ?? throw new Exception("Error in Optimizer. No internal duplication event found.");
            index = events.IndexOf(internalDup);
            newEvents[index] = GetNewEventLength(newEvents[index]);
        }
        var newSignature = new Signature(1, newEvents);
        return currentParams with { Signatures = new Dictionary<string, Signature> { ["CNs"] = newSignature } };
    }
    private CNEventPars GetNewEventLength(CNEventPars oldEvent)
    {
        var newSize = GetNewLength(oldEvent.Size);
        return oldEvent with { Size = newSize };
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

    private double GetMeanSegDistance(Dictionary<string, List<CopyNumber>> simCNPs)
    {
        var obsValues = SummaryFeatures.GetMeanSegLength(ObservedCNPs);
        var simValues = SummaryFeatures.GetMeanSegLength(simCNPs);
        var histMax = Math.Max(obsValues.Max(), simValues.Max());
        var histMin = 0;
        var histBins = 100;
        return CalculateDistance(obsValues, simValues, histBins, histMin, histMax);
    }

    protected double GetPloidyDistance(Dictionary<string, List<CopyNumber>> simCNPs)
    {
        var obsValues = SummaryFeatures.GetPloidy(GenRef, ObservedCNPs, IsFemaleObservedDict, IncludeSexChromosomes);
        var simValues = SummaryFeatures.GetPloidy(GenRef, simCNPs, IsFemaleSimulatedDict, IncludeSexChromosomes);
        var histMax = Math.Max(obsValues.Max(), simValues.Max());
        var histMin = 0;
        var histBins = 100;
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
        //var obsValues = SummaryFeatures.GetBreakpointsDistribution(GenRef, ObservedCNPs, IncludeSexChromosomes, BreakpointsPerChrom, BPBinSize);
        //var simValues = SummaryFeatures.GetBreakpointsDistribution(GenRef, simCNPs, IncludeSexChromosomes, BreakpointsPerChrom, BPBinSize);
        var obsValues = SummaryFeatures.GetBreakpoints(ObservedCNPs, IncludeSexChromosomes);
        var simValues = SummaryFeatures.GetBreakpoints(ObservedCNPs, IncludeSexChromosomes);
        //var histMax = IncludeSexChromosomes ? 24 : 22;
        // Limit the maximum number of breakpoints to 100.
        var histMax = 100;
        var histMin = 0;
        var histBins = histMax;
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