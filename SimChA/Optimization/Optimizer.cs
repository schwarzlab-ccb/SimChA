using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.Simulation;
using SimChA.IO;
using MathNet.Numerics.Statistics;
namespace SimChA.Optimization;

public class Optimizer
{   
    public Dictionary<string, List<CopyNumber>> ObservedCNPs { get; }
    public Dictionary<string, List<CopyNumber>> SimulatedCNPs { get; set;}
    GenRef GenRef { get; }
    public Optimizer(GenRef genRef, List<Sample> observedData)
    {
        GenRef = genRef;
        ObservedCNPs = GetCNPs(GenRef, observedData);
        SimulatedCNPs = new Dictionary<string, List<CopyNumber>>();
    }

    public double Optimize(SimParams simParams, Random rnd, int repeats)
    {
        SimulatedCNPs = GenerateSimulatedCNPs(simParams, rnd, repeats);
        var segDist = GetSegLengthDistance();
        var cpDist = GetChangepointDistance();
        var bpDist = GetBreakpointDistance();
        var majDist = GetMajMinCNDistance(true);
        var minDist = GetMajMinCNDistance(false);
        //Console.WriteLine();
        //Console.WriteLine($"Seg Length WD: {segDist}; Changepoint WD: {cpDist}; BP per chr WD: {bpDist}, maj CNs: {majDist}, min CNs {minDist}");        
        return Math.Sqrt(segDist*segDist + cpDist*cpDist + bpDist*bpDist + majDist*majDist + minDist*minDist);
    }
    public Dictionary<string, List<CopyNumber>> GenerateSimulatedCNPs(SimParams simParams, Random rnd, int repeats)
    {
        if (simParams.Signatures is null || simParams.Signatures.Count == 0)
        {
            throw new Exception("No signatures were provided.");
        }
        Validators.ValidateSignatures(simParams.Signatures);

        var samples = Converters.MakeSamples(rnd, repeats, simParams.EventCount, simParams.EventDist, simParams.Signatures, simParams.Sex, simParams.MCTarget);
        var simulator = new Simulator(rnd, GenRef);
        foreach (var sample in samples)
        {
            simulator.SampleEvents(sample);
        }
        return GetCNPs(GenRef, samples);
    }

    public double GetSegLengthDistance()
    {
        var dataSegList = SummaryFeatures.GetSegLengths(ObservedCNPs);
        var simSegList = SummaryFeatures.GetSegLengths(SimulatedCNPs);
        var histMax = GenRef.ChrLengths["chr1"];
        var histMin = 0;
        var histBins = 200;
        return CalculateDistance(dataSegList, simSegList, histBins, histMin, histMax);
    }

    public double GetChangepointDistance()
    {
        var (obsValues, obsMax) = SummaryFeatures.GetChangepointInfo(ObservedCNPs);
        var (simValues, simMax)  = SummaryFeatures.GetChangepointInfo(SimulatedCNPs);
        var histMax = Math.Max(obsMax, simMax);
        var histMin = 0;
        var histBins = 50;
        return CalculateDistance(obsValues, simValues, histBins, histMin, histMax);

    }

    public double GetBreakpointDistance()
    {
        var (obsValues, obsMax) = SummaryFeatures.GetBreakpointsPerChromosome(ObservedCNPs);
        var (simValues, simMax)  = SummaryFeatures.GetBreakpointsPerChromosome(SimulatedCNPs);
        var histMax = Math.Max(obsMax, simMax);
        var histMin = 0;
        var histBins = 50;
        return CalculateDistance(obsValues, simValues, histBins, histMin, histMax);
    }

    public double GetMajMinCNDistance(bool getMajor)
    {
        var (obsValues, obsMax) = SummaryFeatures.GetMajMinCNs(ObservedCNPs, getMajor);
        var (simValues, simMax) = SummaryFeatures.GetMajMinCNs(SimulatedCNPs, getMajor);
        var histMax = Math.Max(obsMax, simMax);
        var histMin = 0;
        var histBins = 50;
        return CalculateDistance(obsValues, simValues, histBins, histMin, histMax);
    }

    public double CalculateDistance(List<double> data, List<double> sim, int bins, int min, double max)
    {
        var dataHist = new Histogram(data, bins, min, max);
        var simHist  = new Histogram(sim, bins, min, max);
        return StatisticMeasures.WassersteinDistance(dataHist, simHist);
    }


    public Dictionary<string, List<CopyNumber>> GetCNPs(GenRef genRef, List<Sample> samples)
    {
        var cnps = new Dictionary<string, List<CopyNumber>>();
        foreach (var sample in samples)
        {
            foreach (var clone in sample.Clones)
            {
                var cnp = CopyNumbers.CalcCopyNumbers(genRef, sample.Kars[clone.CloneId], sample.Kars[clone.CloneId].SexXX).ToList();
                cnps[sample.SampleId] = cnp;
            }
        }
        return cnps;
    }
}