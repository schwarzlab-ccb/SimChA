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

    private double GetScore(Dictionary<string, List<CopyNumber>> cnps)
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
        var currentParams = GetProposalParams(SimParams);
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
        }
        return bestParams;
    }

    private SimParams GetProposalParams2(SimParams currentParams)
    {

        return currentParams;
    }

    private double GetNewWeight(double oldProb)
    {
        var nTries = 0;
        var newProb = oldProb * (1 + (Rnd.NextDouble() < 0.5 ? -1 : 1) * Rnd.NextDouble() * OptimizationParams.StepFactor);
        while (Math.Abs(newProb - oldProb)/oldProb <= double.Epsilon && nTries < OptimizationParams.MaxTries)
        {
            nTries++;
            newProb = oldProb * (1 + (Rnd.NextDouble() < 0.5 ? -1 : 1) * Rnd.NextDouble() * OptimizationParams.StepFactor);
        }
        return newProb;
    }

    private long GetNewSize(long oldSize)
    {
        var nTries = 0;
        var newSize = (long) (oldSize * (1 + (Rnd.NextDouble() < 0.5 ? -1 : 1) * Rnd.NextDouble() * OptimizationParams.StepFactor));
        while (newSize <= 0 && Math.Abs(newSize-oldSize)/oldSize <= double.Epsilon && nTries < OptimizationParams.MaxTries)
        {
            nTries++;
            newSize = (long) (oldSize * (1 + (Rnd.NextDouble() < 0.5 ? -1 : 1) * Rnd.NextDouble() * OptimizationParams.StepFactor));
        }
        return newSize;
    }
    
    private SimParams GetAllNewParams(SimParams currentParams)
    {
        if (OptimizationParams.ResetSeed)
        {
            currentParams = currentParams with {Seed = -1};
        }
        var events = currentParams.Signatures["CNs"].Events.Where(e => e.Prob > 0).ToList();
        var targetWeight = events.Sum(e => e.Prob);
        var newProbs = new List<double>();
        foreach (var ev in events)
        {
            var weight = GetNewWeight(ev.Prob);
            newProbs.Add(weight);
        }
        var currentTotal = newProbs.Sum();
        newProbs  = newProbs.Select(x => x / currentTotal * targetWeight).ToList();
        var newEvents = events.Select((e, i) => e with { Prob = newProbs[i] }).ToList();

        if (!OptimizationParams.EventWeightsOnly)
        {
            var internalDel = events.Find(e => e.Type == CNEventType.InternalDeletion) ?? throw new Exception("Error in Optimizer. No internal deletion event found.");
            var index = events.IndexOf(internalDel);
            newEvents[index] = newEvents[index] with { Size = GetNewSize(internalDel.Size) };

            var internalDup = events.Find(e => e.Type == CNEventType.InternalDuplication) ?? throw new Exception("Error in Optimizer. No internal duplication event found.");
            index = events.IndexOf(internalDup);
            newEvents[index] = newEvents[index] with { Size = GetNewSize(internalDup.Size) };
        }
        var newSignature = new Signature(1, newEvents);
        return currentParams with { Signatures = new Dictionary<string, Signature> { ["CNs"] = newSignature } };
    }

    private SimParams GetProposalParams(SimParams currentParams)
    {
        var events = currentParams.Signatures["CNs"].Events.Where(e => e.Prob > 0).ToList();
        var totalWeight = events.Sum(e => e.Prob);
        // Reset the seed if selected
        if (OptimizationParams.ResetSeed)
        {
            currentParams = currentParams with {Seed = -1};
        }
        // Choose an event to modify
        // We add an extra two possible indices to account for the InternalDuplication and
        // InternalDeletion Length Parameters
        var nParams = OptimizationParams.EventWeightsOnly ? events.Count : events.Count + 2;
        var index = Rnd.Next(nParams);
        var sign = Rnd.NextDouble() < 0.5 ? -1 : 1;
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
            var newSize = GetNewSize(oldSize);
            var newEvents = new List<CNEventPars>(events);
            newEvents[index] = newEvents[index] with { Size = newSize };
            var newSignature = new Signature(1, newEvents);
            return currentParams with { Signatures = new Dictionary<string, Signature> { ["CNs"] = newSignature } };
        }
        else
        {
            // Modify the relative weight of the event
            var newProb = GetNewWeight(events[index].Prob);
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

    private double GetMeanSegDistance(Dictionary<string, List<CopyNumber>> simCNPs)
    {
        var obsValues = SummaryFeatures.GetMeanSegLength(ObservedCNPs);
        var simValues = SummaryFeatures.GetMeanSegLength(simCNPs);
        var histMax = Math.Max(obsValues.Max(), simValues.Max());
        var histMin = 0;
        var histBins = 100;
        return CalculateDistance(obsValues, simValues, histBins, histMin, histMax);
    }

    private double GetPloidyDistance(Dictionary<string, List<CopyNumber>> simCNPs)
    {
        var (obsValues, obsMax) = SummaryFeatures.GetPloidy(GenRef, ObservedCNPs, IsFemaleObservedDict, IncludeSexChromosomes);
        var (simValues, simMax) = SummaryFeatures.GetPloidy(GenRef, simCNPs, IsFemaleSimulatedDict, IncludeSexChromosomes);
        var histMax = Math.Max(obsMax, simMax);
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
        var obsValues = SummaryFeatures.GetBreakpointsDistribution(GenRef, ObservedCNPs, IncludeSexChromosomes, BreakpointsPerChrom, BPBinSize);
        var simValues = SummaryFeatures.GetBreakpointsDistribution(GenRef, simCNPs, IncludeSexChromosomes, BreakpointsPerChrom, BPBinSize);
        var histMax = IncludeSexChromosomes ? 24 : 22;
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