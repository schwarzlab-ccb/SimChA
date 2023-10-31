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
    public FitnessOptimizer(SimParams simParams, Random rnd, int repeats, 
        GenRef genRef, List<Sample> observedData, List<double> fitnessList) 
        : base(simParams, rnd, repeats, genRef,observedData)
    {
        FitnessList = fitnessList;
        ChromosomeBins = new Dictionary<string, List<long>>();
        SetChromosomeBins();
        ObservedCNPs1MB = GetCNPs(observedData, true);
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
        Validators.ValidateSignatures(SimParams.Signatures);

        var samples = Converters.MakeSamples(Rnd, Repeats, SimParams.EventCount, SimParams.EventDist, SimParams.Signatures, SimParams.Sex, FitnessList);
        var simulator = new Simulator(Rnd, GenRef);
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
        SimulatedCNPs1MB = GetCNPs(samples, true);
        return;
    }

    public double GetFitnessDistance()
    {
        var hdDist = GetHomozygousDeletionDistance();
        var meanCNAcrossGenomeDist = GetMeanCopyNumberAlongGenomeDistance();
        var meanCN = GetMeanCopyNumberDistance();

        return (hdDist + meanCNAcrossGenomeDist + meanCN)/3;
    }

    private void SetChromosomeBins()
    {
        var binSize = 1_000_000;
        foreach (var chrom in GenRef.AllChrs)
        {
            var nFullBins = GenRef.ChrLengths[chrom] / binSize;
            var remainder = GenRef.ChrLengths[chrom] % binSize;
            // Adjusting the first and last bins
            var endBinSize = (long)(0.5 + remainder / 2.0);
            var offset = remainder - 2*endBinSize;
            var binList = new List<long>{0};
            for (int i = 0; i < nFullBins; i++)
            {
                binList.Add(i*binSize+endBinSize);
            }
            binList.Add(binList.Last()+endBinSize+offset-1);
            ChromosomeBins[chrom] = binList;
        }
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


    private Dictionary<string, List<CopyNumber>> GetCNPs(List<Sample> samples, bool binsOf1MB = false)
    {
        var cnps = new Dictionary<string, List<CopyNumber>>();
        foreach (var sample in samples)
        {
            foreach (var clone in sample.Clones)
            {
                var cnp = !binsOf1MB
                    ? CopyNumbers.CalcCopyNumbers(GenRef, sample.Kars[clone.CloneId], sample.Kars[clone.CloneId].SexXX).ToList()
                    : CopyNumbers.CalcCopyNumbers(GenRef, sample.Kars[clone.CloneId], ChromosomeBins, sample.Kars[clone.CloneId].SexXX, true).ToList();
                cnps[sample.SampleId] = cnp;
            }
        }
        return cnps;
    }
}