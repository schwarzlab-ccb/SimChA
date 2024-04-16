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
    protected Dictionary<string, int> ObservedEventCounts { get; }
    protected Dictionary<string, bool> IsFemaleObservedDict {get; set;}
    private bool IncludeSexChromosomes { get; }
    private bool BreakpointsPerChrom {get;}
    private long BPBinSize {get;}
    protected List<Sample> ObservedSamples { get; }
    protected List<Sample> SimulatedSamples { get; set; }
    
    public Optimizer(SimParams simParams, Random rnd, int repeats, GenRef genRef, List<Sample> observedData, bool includeSexChromosomes, SimParams targetParams = null)
    {
        SimParams = simParams;
        Rnd = rnd;
        Repeats = repeats;
        GenRef = genRef;
        ObservedSamples = targetParams == null ? observedData : GenerateSimulatedData(targetParams);
        (ObservedCNPs, ObservedEventCounts) = GetInfo(ObservedSamples);
        IsFemaleObservedDict = ObservedSamples.ToDictionary(s => s.SampleId, s => s.SexXX);
        OptimizationParams = SimParams.OptimizationParams ?? throw new Exception("Error in Optimizer. OptimizationParams not set.");
        IncludeSexChromosomes = includeSexChromosomes;
        BreakpointsPerChrom = OptimizationParams.BreakpointsPerChrom;
        BPBinSize = OptimizationParams.BreakpointsBinSize;
        if (SimParams.Signatures is null || SimParams.Signatures.Count == 0)
        {
            throw new Exception("Error in Optimizer. No signatures were provided.");
        }

        if (OptimizationParams.StepSize <= 0)
        {
            throw new Exception("Error in Optimizer. StepSize must be greater than 0.");
        }
    }

    public virtual SimParams Optimize(FileIO files)
        => OptimizationParams.OptimizationMethod switch
            {
                "MetropolisHastings" => FindBestParams(files),
                "SimulatedAnnealing" => AnnealingBestParams(files),
                "AdaptiveSimulatedAnnealing" => AnnealingBestParams(files),
                "StepSizeDecay" => DecayBestParams(files),
                _ => throw new Exception("Error in Optimizer. Optimization method not recognized."),
            };

    public double GetABCDistance()
    {
        var samples = GenerateSimulatedData(SimParams);
        var (simCNPs, eventCount) = GetInfo(samples);
        return GetScore(samples);
    }

    public double GetScore(List<Sample> samples)
    {
        var (cnps, eventCounts) = GetInfo(samples);
        var isFemaleDict = samples.ToDictionary(s => s.SampleId, s => s.SexXX);
        var totalDist = new List<double>();
        if (OptimizationParams.UseSegLength)
        {
            var segDist = GetSegLengthDistance(cnps);
            totalDist.Add(segDist*segDist);
        }
        if (OptimizationParams.UsePloidy)
        {
            var ploidyDist = GetPloidyDistance(cnps, isFemaleDict);
            totalDist.Add(ploidyDist*ploidyDist);
        }
        if (OptimizationParams.UseBreakpoints)
        {
            var bpDist = GetBreakpointDistance(cnps, eventCounts);
            totalDist.Add(bpDist*bpDist);
        }
        //var copyNumberMatrix = SummaryFeatures.GetChrCopyNumberMatrix(GenRef.AllChrs, cnps);
        //var mkv = SummaryFeatures.GetMKV(copyNumberMatrix);
        //var aneuploidy = SummaryFeatures.GetAverageAneuploidy(copyNumberMatrix);
        return totalDist.Sum();
    }

    private double GetAcceptanceProbability(double scoreA, double scoreB, double temperature)
        => Math.Min(0, -OptimizationParams.AcceptanceFactor*(scoreB - scoreA)/(scoreA*temperature));

    protected double GetAcceptanceProbability(double scoreA, double scoreB)
        => GetAcceptanceProbability(scoreA, scoreB, 1.0);

    private SimParams FindBestParams(FileIO files)
    {
        var currentParams = SimParams;
        var currentSamples = GenerateSimulatedData(currentParams);
        var currentScore = GetScore(currentSamples);
        var bestParams = currentParams;
        var bestScore = currentScore;
        var counter = 0;
        var stepSize = OptimizationParams.StepSize;
        for (int i = 0; i < OptimizationParams.NumSamplesTotal; i++)
        {
            //Console.WriteLine($"Iteration {i+1} of {OptimizationParams.NumSamplesTotal}");
            var proposedParams = GetProposalParams(currentParams, stepSize);
            var proposedSamples = GenerateSimulatedData(proposedParams);
            var proposedScore = GetScore(proposedSamples);
            var prob = GetAcceptanceProbability(currentScore, proposedScore);
            if (prob >= Math.Log(Rnd.NextDouble()))
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

    private SimParams DecayBestParams(FileIO files)
    {
        var currentParams = SimParams;
        var currentSamples = GenerateSimulatedData(currentParams);
        var currentScore = GetScore(currentSamples);
        var bestParams = currentParams;
        var bestScore = currentScore;
        var counter = 0;
        var stepSize = OptimizationParams.StepSize;
        var decay = OptimizationParams.CoolingRate;

        for (int i = 0; i < OptimizationParams.NumSamplesTotal; i++)
        {
            Console.WriteLine($"Iteration {i+1} of {OptimizationParams.NumSamplesTotal}");
            var proposedParams = GetProposalParams(currentParams, stepSize);
            var proposedSamples = GenerateSimulatedData(proposedParams);
            var proposedScore = GetScore(proposedSamples);
            var prob = GetAcceptanceProbability(currentScore, proposedScore);
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
            // Apply decay to step size
            switch (OptimizationParams.StepSizeDecayType)
            {
                case "Exponential":
                    stepSize *= decay;
                    break;
                case "Linear":
                    stepSize = Math.Max(OptimizationParams.MinStepSize, stepSize - decay);
                    break;
                case "Inverse":
                    stepSize = OptimizationParams.StepSize / (1 + decay * i);
                    break;
                default:
                    throw new Exception("Error in Optimizer. StepSizeDecayType not recognized.");
            }

        }
        return bestParams;
    }

    private SimParams AnnealingBestParams(FileIO files)
    {
        var currentParams = SimParams;
        var currentSamples = GenerateSimulatedData(currentParams);
        var currentScore = GetScore(currentSamples);
        var bestParams = currentParams;
        var bestScore = currentScore;
        var counter = 0;
        // Simulated Annealing parameters
        var temperature = OptimizationParams.StartTemp;
        var coolingRate = OptimizationParams.CoolingRate;
        // Adaptive Simulated Annealing parameters
        var nSuccesses = 0;
        var nFailures = 0;
        var stepSize = OptimizationParams.StepSize;

        for (int i = 0; i < OptimizationParams.NumSamplesTotal; i++)
        {
            // Simulated Annealing updates temperature at the beginning of each iteration
            if (OptimizationParams.OptimizationMethod == "SimulatedAnnealing")
            {
                temperature = OptimizationParams.StartTemp * (1.0 - i/(double)OptimizationParams.NumSamplesTotal);
            }
            Console.WriteLine($"Iteration {i+1} of {OptimizationParams.NumSamplesTotal}");
            var proposedParams = GetProposalParams(currentParams, stepSize);
            var proposedSamples = GenerateSimulatedData(proposedParams);
            var proposedScore = GetScore(proposedSamples);
            var prob = GetAcceptanceProbability(currentScore, proposedScore, temperature);
            if (Rnd.NextDouble() < prob)
            {
                currentParams = proposedParams;
                currentScore = proposedScore;
                nSuccesses++;
            }
            else 
            {
                nFailures++;
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
            // Adaptive Simulated Annealing updates temperature at the end of each iteration
            // Adaptive adjustments

            if (OptimizationParams.OptimizationMethod == "AdaptiveSimulatedAnnealing")
            {
                if (nSuccesses > nFailures)
                {
                    stepSize *= 1.05; // Increase step size to explore more aggressively
                    temperature /= coolingRate; // Cool slower if we're having success
                }
                else
                {
                    // Decrease step size to refine the search
                    stepSize = Math.Max(OptimizationParams.MinStepSize, stepSize*0.95); 
                    temperature *= coolingRate; // Cool faster if we're stuck
                }
                
                // Reset success and failure counters periodically
                if ((i + 1) % 100 == 0)
                {
                    nSuccesses = 0;
                    nFailures = 0;
                }
            }
        }
        return bestParams;
    }

    private SimParams GetProposalParams(SimParams currentParams, double stepSize)
        => OptimizationParams.ParamVariationMode switch
            {
                0 => GetAllNewParams(currentParams, stepSize),
                _ => GetNNewParams(currentParams, stepSize)
            };

    private SimParams GetNNewParams(SimParams currentParams, double stepSize)
    {
        if (OptimizationParams.ResetSeed)
        {
            currentParams = currentParams with { Seed = -1 };
        }
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
                newEvents[eventIndex] = GetNewEventLength(newEvents[eventIndex], stepSize);
            }
            else
            {   
                newProbs[index] = GetNewWeight(events[index].Prob, stepSize);
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

    protected double GetNewWeight(double oldWeight, double stepSize, double minimum = 0.0)
    {
        var nTries = 0;
        double sign = Rnd.NextDouble() < 0.5 ? -1 : 1;
        var newWeight = oldWeight * (1 + sign * Rnd.NextDouble() * stepSize);
        while (newWeight <= minimum && nTries < OptimizationParams.MaxTries)
        {
            nTries++;
            sign = Rnd.NextDouble() < 0.5 ? -1 : 1;
            newWeight = oldWeight * (1 + sign * Rnd.NextDouble() * stepSize);
        }
        // If we don't find a valid value, return the old one 
        if (newWeight < minimum)
        {
            return oldWeight;
        }
        return newWeight;
    }

    private long GetNewLength()
        => Rnd.Next(OptimizationParams.MinLength, OptimizationParams.MaxLength);

    // Lengths have to be greater than or equal to 1
    private long GetNewLength(long oldValue, double stepSize)
     => OptimizationParams.BoundedLengths ? GetNewLength() : (long)GetNewWeight(oldValue, stepSize, 1.0);
    
    private SimParams GetAllNewParams(SimParams currentParams, double stepSize)
    {
        if (OptimizationParams.ResetSeed)
        {
            currentParams = currentParams with { Seed = -1 };
        }
        var events = currentParams.Signatures["CNs"].Events.Where(e => e.Prob > 0).ToList();
        var targetWeight = events.Sum(e => e.Prob);
        var newProbs = new List<double>();
        foreach (var ev in events)
        {
            var weight = GetNewWeight(ev.Prob, stepSize);
            newProbs.Add(weight);
        }
        var currentTotal = newProbs.Sum();
        newProbs  = newProbs.Select(x => targetWeight * x / currentTotal).ToList();
        var newEvents = events.Select((e, i) => e with { Prob = newProbs[i] }).ToList();

        if (!OptimizationParams.EventWeightsOnly)
        {
            var internalDel = events.Find(e => e.Type == CNEventType.InternalDeletion) ?? throw new Exception("Error in Optimizer. No internal deletion event found.");
            var index = events.IndexOf(internalDel);
            newEvents[index] = GetNewEventLength(newEvents[index], stepSize);

            var internalDup = events.Find(e => e.Type == CNEventType.InternalDuplication) ?? throw new Exception("Error in Optimizer. No internal duplication event found.");
            index = events.IndexOf(internalDup);
            newEvents[index] = GetNewEventLength(newEvents[index], stepSize);
        }
        var newSignature = new Signature(1, newEvents);
        return currentParams with { Signatures = new Dictionary<string, Signature> { ["CNs"] = newSignature } };
    }
    private CNEventPars GetNewEventLength(CNEventPars oldEvent, double stepSize)
    {
        var newSize = GetNewLength(oldEvent.Size, stepSize);
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
    protected double GetSegLengthDistance(Dictionary<string, List<CopyNumber>> simCNPs)
        => OptimizationParams.SegLengthType switch
            {
                "Stratified" => GetStratifiedSegLengthDistance(simCNPs),
                "All" => GetAllSegLengths(simCNPs),
                "Mean" => GetMeanSegDistance(simCNPs),
                _ => throw new Exception("Error in Optimizer. SegLengthType not recognized."),
            };

    private double GetAllSegLengths(Dictionary<string, List<CopyNumber>> simCNPs)
    {
        long cutoff = -1;
        var histMax = GenRef.ChrLengths["chr1"] + 1.0;
        var histMin = 0.0;
        var histBins = 101;//(int) (histMax / 1_000_000);
        if (OptimizationParams.SegLengthCutoff > 0)
        {
            cutoff = OptimizationParams.SegLengthCutoff;
            histMax = (int)cutoff;
            histBins = 100;
        }
        var obsValues = SummaryFeatures.GetSegLengths(ObservedCNPs, cutoff);
        var simValues = SummaryFeatures.GetSegLengths(simCNPs, cutoff);
        if (OptimizationParams.LogTransformSegLength)
        {
            obsValues = obsValues.Select(x => Math.Log(x)).ToList();
            simValues = simValues.Select(x => Math.Log(x)).ToList();
            histMax = Math.Log(histMax);
            histMin = Math.Log(1/1000000.0);
        }
        return CalculateDistance(obsValues, simValues, histBins, histMin, histMax);
    }

    private double GetStratifiedSegLengthDistance(Dictionary<string, List<CopyNumber>> simCNPs)
    {
        var histMax = (double)GenRef.ChrLengths["chr1"];
        var histMin = 0.0;
        var histBins = 101;
        var weighted = OptimizationParams.SegmentCountWeighted;
        var obsValues = SummaryFeatures.GetStratifiedSegLengths(ObservedCNPs, weighted);
        var simValues = SummaryFeatures.GetStratifiedSegLengths(simCNPs, weighted);
        if (OptimizationParams.LogTransformSegLength)
        {
            for (int i = 0; i < obsValues.Count; i++)
            {
                obsValues[i] = (obsValues[i].weight, obsValues[i].segs.Select(x => Math.Log(x)).ToList());
                simValues[i] = (simValues[i].weight, simValues[i].segs.Select(x => Math.Log(x)).ToList());
            }
            histMax = Math.Log(histMax) + 1.0;
            histMin = Math.Log(1/1000000.0);
        }
        var totalDist = new List<double>();
        for (int i = 0; i < obsValues.Count; i++)
        {
            var obs = obsValues[i].segs;
            var sim = simValues[i].segs;
            totalDist.Add(CalculateDistance(obs, sim, histBins, histMin, histMax) * simValues[i].weight);
        }
        return totalDist.Sum();
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

    protected double GetPloidyDistance(Dictionary<string, List<CopyNumber>> simCNPs, Dictionary<string, bool> simFemaleDict)
    {
        double cutoff = -1.0;
        if (OptimizationParams.UsePloidy)
        {
            cutoff = OptimizationParams.PloidyCutoff;
        }
        var obsValues = SummaryFeatures.GetPloidy(GenRef, ObservedCNPs, IsFemaleObservedDict, IncludeSexChromosomes, cutoff);
        var simValues = SummaryFeatures.GetPloidy(GenRef, simCNPs, simFemaleDict, IncludeSexChromosomes, cutoff);
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
    protected double GetBreakpointDistance(Dictionary<string, List<CopyNumber>> simCNPs, Dictionary<string, int> simEventCounts)
    {
        //var obsValues = SummaryFeatures.GetBreakpointsDistribution(GenRef, ObservedCNPs, IncludeSexChromosomes, BreakpointsPerChrom, BPBinSize);
        //var simValues = SummaryFeatures.GetBreakpointsDistribution(GenRef, simCNPs, IncludeSexChromosomes, BreakpointsPerChrom, BPBinSize);
        var obsValues = SummaryFeatures.GetBreakpoints(ObservedCNPs, ObservedEventCounts, IncludeSexChromosomes);
        var simValues = SummaryFeatures.GetBreakpoints(simCNPs, simEventCounts, IncludeSexChromosomes);
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
    protected static double CalculateDistance(List<double> data, List<double> sim, int bins, double min, double max)
    {
        var dataHist = new Histogram(data, bins, min, max);
        var simHist  = new Histogram(sim, bins, min, max);
        return StatisticMeasures.WassersteinDistance(dataHist, simHist);
    }
    protected (Dictionary<string, List<CopyNumber>> cnps, Dictionary<string, int> eventCounts) GetInfo(List<Sample> samples)
    {
        var cnps = new Dictionary<string, List<CopyNumber>>();
        var eventCounts = new Dictionary<string, int>();
        foreach (var sample in samples)
        {
            foreach (var clone in sample.Clones)
            {
                var cnp = CopyNumbers.CalcCopyNumbers(GenRef, sample.Kars[clone.CloneId], sample.Kars[clone.CloneId].SexXX).ToList();
                cnps[sample.SampleId] = cnp;
                eventCounts[sample.SampleId] = sample.EventDescs[clone.CloneId].Count;
            }
        }
        return (cnps, eventCounts);
    }
}