using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.Simulation;
using SimChA.IO;
using MathNet.Numerics.Statistics;
namespace SimChA.Optimization;

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

    public virtual double Optimize()
    {
        GenerateSimulatedCNPs();
        return GetEventDistance();
    }

    private double GetEventDistance()
    {
        var segDist = GetSegLengthDistance();
        var cpDist = GetChangepointDistance();
        var bpDist = GetBreakpointDistance();
        var majDist = GetMajMinCNDistance(true);
        var minDist = GetMajMinCNDistance(false);      
        //return Math.Sqrt(segDist*segDist + cpDist*cpDist + bpDist*bpDist + majDist*majDist + minDist*minDist);
        return (segDist + cpDist + bpDist + majDist + minDist)/5;
    }
    private List<Sample> GenerateSimulatedData()
    {
        if (SimParams.Signatures is null || SimParams.Signatures.Count == 0)
        {
            throw new Exception("No signatures were provided.");
        }
        Validators.ValidateSignatures(SimParams.Signatures);

        var samples = Converters.MakeSamples(Rnd, Repeats, SimParams.EventCount, SimParams.EventDist, SimParams.Signatures, SimParams.Sex, SimParams.MCTarget);
        var simulator = new Simulator(Rnd, GenRef);
        foreach (var sample in samples)
        {
            simulator.SampleEvents(sample);
        }
        return samples;
    }
    private void GenerateSimulatedCNPs()
    {
        var samples = GenerateSimulatedData();
        SimulatedCNPs = GetCNPs(samples);
        return;
    }
    private double GetSegLengthDistance()
    {
        var (obsValues, obsMax) = SummaryFeatures.GetSegLengths(ObservedCNPs);
        var (simValues, simMax) = SummaryFeatures.GetSegLengths(SimulatedCNPs);
        var histMax = Math.Max(obsMax, simMax);
        var histMin = 0;
        var histBins = 200;
        return CalculateDistance(obsValues, simValues, histBins, histMin, histMax);
    }

    private double GetChangepointDistance()
    {
        var (obsValues, obsMax) = SummaryFeatures.GetChangepointInfo(ObservedCNPs);
        var (simValues, simMax)  = SummaryFeatures.GetChangepointInfo(SimulatedCNPs);
        var histMax = Math.Max(obsMax, simMax);
        var histMin = 0;
        var histBins = 50;
        return CalculateDistance(obsValues, simValues, histBins, histMin, histMax);
    }
    private double GetBreakpointDistance()
    {
        var (obsValues, obsMax) = SummaryFeatures.GetBreakpointsPerChromosome(ObservedCNPs);
        var (simValues, simMax)  = SummaryFeatures.GetBreakpointsPerChromosome(SimulatedCNPs);
        var histMax = Math.Max(obsMax, simMax);
        var histMin = 0;
        var histBins = 50;
        return CalculateDistance(obsValues, simValues, histBins, histMin, histMax);
    }
    private double GetMajMinCNDistance(bool getMajor)
    {
        var (obsValues, obsMax) = SummaryFeatures.GetMajMinCNs(ObservedCNPs, getMajor);
        var (simValues, simMax) = SummaryFeatures.GetMajMinCNs(SimulatedCNPs, getMajor);
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
    private Dictionary<string, List<CopyNumber>> GetCNPs(List<Sample> samples)
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