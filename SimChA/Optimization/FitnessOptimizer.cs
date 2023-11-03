using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.Simulation;
using SimChA.IO;
using MathNet.Numerics.Statistics;
namespace SimChA.Optimization;

public class FitnessOptimizer : Optimizer
{   
    private Dictionary<string, List<CopyNumber>> ObservedCNPs { get; }
    private Dictionary<string, List<CopyNumber>> ObservedCNPs1MB { get; }
    private Dictionary<string, List<CopyNumber>> SimulatedCNPs { get; set;}
    private Dictionary<string, List<CopyNumber>> SimulatedCNPs1MB { get; set;}
    protected readonly Dictionary<string, List<long>> ChromosomeBins;
    private readonly List<double> FitnessList;
    private Binner Binner;
    public FitnessOptimizer(SimParams simParams, Random rnd, int repeats, 
        GenRef genRef, List<Sample> observedData, string binnedSamples, List<double> fitnessList) 
        : base(simParams, rnd, repeats, genRef, observedData)
    {
        FitnessList = fitnessList;
        Binner = new Binner(GenRef);
        ObservedCNPs1MB = FileIO.ReadProfiles(binnedSamples);
        ObservedCNPs    = GetCNPs(observedData);
        SimulatedCNPs = new Dictionary<string, List<CopyNumber>>();
        SimulatedCNPs1MB = new Dictionary<string, List<CopyNumber>>();
    }

    public override double Optimize()
    {
        GenerateSimulatedCNPs();
        return GetFitnessDistance();
    }

    private List<Sample> GenerateSimulatedData()
    {
        if (SimParams.Signatures is null || SimParams.Signatures.Count == 0)
        {
            throw new Exception("No signatures were provided.");
        }
        if (SimParams.MCParams is null)
        {
            throw new Exception("No MC parameters were provided.");
        }
        Validators.ValidateSignatures(SimParams.Signatures);

        var samples = Converters.MakeSamples(Rnd, Repeats, SimParams.EventCount, SimParams.EventDist, SimParams.Signatures, SimParams.Sex, FitnessList);
        var simulator = new MCSimulator(Rnd, GenRef, SimParams.Fitness, SimParams.MCParams);
        foreach (var sample in samples)
        {
            simulator.SampleEvents(sample);
        }
        return samples;
    }
    public virtual void GenerateSimulatedCNPs()
    {
        var samples = GenerateSimulatedData();
        SimulatedCNPs = GetCNPs(samples);
        SimulatedCNPs1MB = Binner.GetBinnedCNProfiles(samples);
        return;
    }

    public double GetFitnessDistance()
    {
        var hdDist = GetHomozygousDeletionDistance();
        var meanCNAcrossGenomeDist = GetMeanCopyNumberAlongGenomeDistance();
        var meanCN = GetMeanCopyNumberDistance();

        return (hdDist + meanCNAcrossGenomeDist + meanCN)/3;
    }

    private double GetMeanCopyNumberAlongGenomeDistance()
    {
        // TODO: Do we worry about the slightly smaller bins from the fact that the genome 
        // is not completely divisible by 1MB?
        // What about the segments that are partially in the binned region?
        var obsCounts = SummaryFeatures.GetMeanCopyNumberAlongGenome(ObservedCNPs1MB);
        var simCounts = SummaryFeatures.GetMeanCopyNumberAlongGenome(SimulatedCNPs1MB);
        
        return StatisticMeasures.WassersteinDistance(obsCounts, simCounts);
    }

    private double GetMeanCopyNumberDistance()
    {
        var obsValues = ObservedCNPs.SelectMany(cnp => cnp.Value)
                        .Where(cn => cn.CNH1 + cn.CNH2 >= 0)
                        .Select(cn => (double)cn.CNH1 + cn.CNH2).Average();
        var simValues = SimulatedCNPs.SelectMany(cnp => cnp.Value)
                        .Select(cn => (double)cn.CNH1 + cn.CNH2).Average();
        return Math.Abs(obsValues - simValues);
    }

    private double GetHomozygousDeletionDistance()
    {
        var (obsValues, obsMax) = SummaryFeatures.GetHomozygousDeletionFraction(ObservedCNPs, GenRef.AutosomeLen);
        var (simValues, simMax)  = SummaryFeatures.GetHomozygousDeletionFraction(SimulatedCNPs, GenRef.AutosomeLen);
        var histMax = Math.Max(obsMax, simMax);
        var histMin = 0;
        var histBins = 50;
        return CalculateDistance(obsValues, simValues, histBins, histMin, histMax);
    }
}