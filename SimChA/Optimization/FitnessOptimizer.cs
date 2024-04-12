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
    private Dictionary<string, List<CopyNumber>> ObservedCNPs1MB { get; }

    private readonly List<(double, double, double)> FitnessComponents;
    private readonly Binner Binner;

    public FitnessOptimizer(SimParams simParams, Random rnd, int repeats, 
        GenRef genRef, List<Sample> observedData, bool includeSexChromosomes, string binnedSamples, List<(double, double, double)> fitnessComponents) 
        : base(simParams, rnd, repeats, genRef, observedData, includeSexChromosomes)
    {
        FitnessComponents = fitnessComponents;
        Binner = new Binner(GenRef);
        ObservedCNPs1MB = FileIO.ReadProfiles(binnedSamples);
    }

    public override SimParams Optimize(FileIO files)
        => FindBestParams(files);

    /*private double GetScore(List<Sample> samples)
    {
        var cnps = GetCNPs(samples);
        var binnedCNPs = Binner.GetBinnedCNProfiles(samples);
        IsFemaleSimulatedDict = samples.ToDictionary(s => s.SampleId, s => s.SexXX);
        var distance = GetFitnessDistance(cnps, binnedCNPs);
        return distance;
    }*/

    private SimParams FindBestParams(FileIO files)
    {
        var currentParams = SimParams;
        var currentSamples = GenerateSimulatedData(currentParams);
        var currentScore = GetScore(currentSamples);
        var bestParams = currentParams;
        var bestScore = currentScore;
        var counter = 0;
        var stepSize = OptimizationParams.StepSize;
        if (stepSize <= 0)
        {
            throw new Exception("Error in FitnessOptimizer. Step size must be greater than 0.");
        }
        for (int i = 0; i < OptimizationParams.NumSamplesTotal; i++)
        {
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
                //files.WriteSimParams(bestParams, $"best_params_{counter}.json");
                files.WriteSimParams(bestParams);
                counter++;
            }
             if (OptimizationParams.WriteIntermediate && i % OptimizationParams.WriteFrequency == 0)
            {
                files.WriteSimParams(currentParams, $"params_{i}.json");
            }
        }
        return bestParams;
    }

    private SimParams GetAllNewParams(SimParams currentParams, double stepSize)
    {
        var fitnessList = currentParams.Fitness.ParamsList();
        var newParams = new List<double>();
        foreach (var oldValue in fitnessList)
        {
            newParams.Add(GetNewWeight(oldValue, stepSize));
        }
        var newFitness = new FitnessParams(newParams[0], newParams[1], newParams[2], newParams[3]);
        return currentParams with { Fitness = newFitness };
    }

    private SimParams GetOneNewParam(SimParams currentParams, double stepSize)
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
        double newWeight = oldWeight * (1 + sign * Rnd.NextDouble() * stepSize);
        int nTries = 0;
        while (Math.Abs(newWeight - oldWeight)/oldWeight <= double.Epsilon && nTries < 10)
        {
            nTries++;
            sign = Rnd.NextDouble() < 0.5 ? -1 : 1;
            newWeight = oldWeight * (1 + sign * Rnd.NextDouble() * OptimizationParams.StepSize);
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

    private SimParams GetProposalParams(SimParams currentParams, double stepSize)
        => OptimizationParams.ParamVariationMode switch
            {
                0 => GetAllNewParams(currentParams, stepSize),
                1 => GetOneNewParam(currentParams, stepSize),
                _ => throw new Exception("Error in FitnessOptimizer. ParamVariationMode is not valid.")
            };
    private List<double> GetFitnessList(FitnessParams fParams)
    {
        var fitnessList = new List<double>();
        foreach (var component in FitnessComponents)
        {
            double stressTerm = component.Item1 * fParams.Stress;
            double tsgogTerm  = component.Item2 * fParams.TsgOg;
            double essTerm    = component.Item3 * fParams.Essentiality;
            fitnessList.Add(1.0 + (stressTerm + tsgogTerm + essTerm)*fParams.TotalStrength);
        }
        return fitnessList;
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
        var fitnessList = GetFitnessList(currentParams.Fitness);
        var samples = Converters.MakeSamples(Rnd, Repeats, currentParams.EventCount, 
            currentParams.EventDist, currentParams.Signatures, currentParams.Sex, fitnessList);
        var simulator = new MCSimulator(Rnd, GenRef, currentParams.Fitness, currentParams.MCParams);
        foreach (var sample in samples)
        {
            simulator.SampleEvents(sample);
        }
        return samples;
    }

    public double GetFitnessDistance(Dictionary<string, List<CopyNumber>> cnps, Dictionary<string, List<CopyNumber>> binnedCNPs)
    {
        var totalDist = new List<double>();
        if (OptimizationParams.UseSegLength)
        {
            var segDist = GetSegLengthDistance(cnps);
            totalDist.Add(segDist*segDist);
        }
        if (OptimizationParams.UseCNAlongGenome)
        {
            var dist = GetMeanCNAlongGenomeDistance(binnedCNPs);
            totalDist.Add(dist*dist);
        }
        if (OptimizationParams.UseBreakpoints)
        {
            var bpDist = GetBreakpointDistance(cnps, new Dictionary<string, int>());
            totalDist.Add(bpDist*bpDist);
        }
        if (OptimizationParams.UseHomozygousDeletion)
        {
            var hdDist = GetHomozygousDeletionDistance(cnps);
            totalDist.Add(hdDist*hdDist);
        }
        if (OptimizationParams.UsePloidy)
        {
            var ploidyDist = GetPloidyDistance(cnps, new Dictionary<string, bool> ());
            totalDist.Add(ploidyDist*ploidyDist);
        }

        return totalDist.Sum();
    }

    private double GetMeanCNAlongGenomeDistance(Dictionary<string, List<CopyNumber>> binnedCNPs)
    {
        // TODO: Do we worry about the slightly smaller bins from the fact that the genome 
        // is not completely divisible by 1MB?
        // What about the segments that are partially in the binned region?
        var obsCounts = SummaryFeatures.GetMeanCNAlongGenome(ObservedCNPs1MB);
        var simCounts = SummaryFeatures.GetMeanCNAlongGenome(binnedCNPs);
        
        return StatisticMeasures.WassersteinDistance(obsCounts, simCounts);
    }

    private double GetHomozygousDeletionDistance(Dictionary<string, List<CopyNumber>> cnps)
    {
        var obsValues = SummaryFeatures.GetHomozygousDeletionFraction(GenRef, ObservedCNPs);
        var simValues  = SummaryFeatures.GetHomozygousDeletionFraction(GenRef, cnps);
        var histMax = Math.Max(obsValues.Max(), simValues.Max());
        var histMin = 0;
        var histBins = 50;
        return CalculateDistance(obsValues, simValues, histBins, histMin, histMax);
    }
}