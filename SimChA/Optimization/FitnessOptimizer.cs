using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.Simulation;
using SimChA.IO;
using MathNet.Numerics.Statistics;
namespace SimChA.Optimization;

using System.CodeDom.Compiler;
using System.Reflection;

public class FitnessOptimizer : Optimizer
{   
    private Dictionary<string, List<CopyNumber>> ObservedCNPs { get; }
    private Dictionary<string, List<CopyNumber>> ObservedCNPs1MB { get; }

    private Dictionary<string, bool> IsFemaleObservedDict {get; set;} 
    private readonly List<double> FitnessList;
    private readonly Binner Binner;
    public FitnessOptimizer(SimParams simParams, Random rnd, int repeats, 
        GenRef genRef, List<Sample> observedData, string binnedSamples, List<double> fitnessList) 
        : base(simParams, rnd, repeats, genRef, observedData)
    {
        FitnessList = fitnessList;
        Binner = new Binner(GenRef);

        ObservedCNPs1MB = FileIO.ReadProfiles(binnedSamples);
        ObservedCNPs    = GetCNPs(observedData);
        IsFemaleObservedDict = observedData.ToDictionary(s => s.SampleId, s => s.SexXX);
    }

    public override SimParams Optimize(FileIO files)
        => FindBestParams(files, 5000, 0.01); // 1000 samples, 1% step size

    private double GetScore(List<Sample> samples)
    {
        var cnps = GetCNPs(samples);
        var binnedCNPs = Binner.GetBinnedCNProfiles(samples);
        var isFemaleDict = samples.ToDictionary(s => s.SampleId, s => s.SexXX);
        var distance = GetFitnessDistance(cnps, binnedCNPs, isFemaleDict);
        return distance;
    }
    private SimParams FindBestParams(FileIO files, int numSamples, double stepFactor)
    {
        var currentParams = SimParams;
        var currentSamples = GenerateSimulatedData(currentParams);
        var currentScore = GetScore(currentSamples);
        for (int i = 0; i < numSamples; i++)
        {
            var proposedParams = GetProposalParams(currentParams, stepFactor);
            var proposedSamples = GenerateSimulatedData(proposedParams);
            var proposedScore = GetScore(proposedSamples);
            var delta = proposedScore - currentScore;
            var prob = Math.Min(1, Math.Exp(-delta));
            if (Rnd.NextDouble() < prob)
            {
                currentParams = proposedParams;
                currentScore = proposedScore;
                files.WriteSimParams(currentParams);
            }
        }
        return currentParams;
    }
    private SimParams GetProposalParams(SimParams currentParams, double stepFactor)
    {
        if (currentParams.MCParams is null)
        {
            throw new Exception("Error in FitnessOptimizer. No MC parameters were provided.");
        }
        int index = Rnd.Next(4);
        // Modify the relative weight of the fitness parameter
        var sign = Rnd.NextDouble() < 0.5 ? -1 : 1;
        double oldWeight = 1;
        var oldStress = currentParams.Fitness.Stress;
        var oldTsgOg = currentParams.Fitness.TsgOg;
        var oldEssentiality = currentParams.Fitness.Essentiality;
        var oldStrength = currentParams.Fitness.TotalStrength;
        switch (index)
        {
            case 0:
                oldWeight = oldStress;
                break;
            case 1:
                oldWeight *= oldTsgOg;
                break;
            case 2:
                oldWeight *= oldEssentiality;
                break;
            default:
                oldWeight *= oldStrength;
                break;
        }
        double newWeight = oldWeight * (1 + sign * Rnd.NextDouble()* stepFactor);
        int nTries = 0;
        while (Math.Abs(newWeight - oldWeight)/oldWeight <= double.Epsilon && nTries < 10)
        {
            nTries++;
            sign = Rnd.NextDouble() < 0.5 ? -1 : 1;
            newWeight = oldWeight * (1 + sign * Rnd.NextDouble() * stepFactor);
        }
        // Stress, TSG/OG, and Essentiality must sum to 1, so if one changed, then
        // the others must change as well.
        var newFactor = (1 - newWeight)/(1-oldWeight);
        return index switch
        {
            0 => currentParams with { Fitness = new FitnessParams(newWeight, oldTsgOg * newFactor, oldEssentiality * newFactor, oldStrength) },
            1 => currentParams with { Fitness = new FitnessParams(oldStress * newFactor, newWeight, oldEssentiality * newFactor, oldStrength) },
            2 => currentParams with { Fitness = new FitnessParams(oldStress * newFactor, oldTsgOg * newFactor, newWeight, oldStrength) },
            _ => currentParams with { Fitness = new FitnessParams(oldStress, oldTsgOg, oldEssentiality, newWeight) },
        };
    }

    private List<Sample> GenerateSimulatedData(SimParams currentParams)
    {
        if (currentParams.Signatures is null || currentParams.Signatures.Count == 0)
        {
            throw new Exception("No signatures were provided.");
        }
        if (currentParams.MCParams is null)
        {
            throw new Exception("No MC parameters were provided.");
        }
        Validators.ValidateSignatures(currentParams.Signatures);

        var samples = Converters.MakeSamples(Rnd, Repeats, currentParams.EventCount, 
            currentParams.EventDist, currentParams.Signatures, currentParams.Sex, FitnessList);
        var simulator = new MCSimulator(Rnd, GenRef, currentParams.Fitness, currentParams.MCParams);
        foreach (var sample in samples)
        {
            simulator.SampleEvents(sample);
        }
        return samples;
    }

    public double GetFitnessDistance(Dictionary<string, List<CopyNumber>> cnps, Dictionary<string, List<CopyNumber>> binnedCNPs, Dictionary<string, bool> isFemaleDict)
    {
        var hdDist = GetHomozygousDeletionDistance(cnps);
        var meanCNAcrossGenomeDist = GetMeanCopyNumberAlongGenomeDistance(binnedCNPs);
        var meanCN = GetPloidyDistance(cnps, isFemaleDict);

        return (hdDist + meanCNAcrossGenomeDist + meanCN)/3;
    }

    private double GetMeanCopyNumberAlongGenomeDistance(Dictionary<string, List<CopyNumber>> binnedCNPs)
    {
        // TODO: Do we worry about the slightly smaller bins from the fact that the genome 
        // is not completely divisible by 1MB?
        // What about the segments that are partially in the binned region?
        var obsCounts = SummaryFeatures.GetMeanCopyNumberAlongGenome(ObservedCNPs1MB);
        var simCounts = SummaryFeatures.GetMeanCopyNumberAlongGenome(binnedCNPs);
        
        return StatisticMeasures.WassersteinDistance(obsCounts, simCounts);
    }

    private double GetPloidyDistance(Dictionary<string, List<CopyNumber>> cnps, Dictionary<string, bool> isFemaleDict)
    {
        var (obsValues, obsMax) = SummaryFeatures.GetPloidy(GenRef, ObservedCNPs, IsFemaleObservedDict);
        var (simValues, simMax) = SummaryFeatures.GetPloidy(GenRef, cnps, isFemaleDict);
        var histMax = Math.Max(obsMax, simMax);
        var histMin = 0;
        var histBins = 50;
        return CalculateDistance(obsValues, simValues, histBins, histMin, histMax);
    }

    private double GetHomozygousDeletionDistance(Dictionary<string, List<CopyNumber>> cnps)
    {
        var (obsValues, obsMax) = SummaryFeatures.GetHomozygousDeletionFraction(ObservedCNPs, GenRef.AutosomeLen);
        var (simValues, simMax)  = SummaryFeatures.GetHomozygousDeletionFraction(cnps, GenRef.AutosomeLen);
        var histMax = Math.Max(obsMax, simMax);
        var histMin = 0;
        var histBins = 50;
        return CalculateDistance(obsValues, simValues, histBins, histMin, histMax);
    }
}