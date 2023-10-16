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
        //ObservedSegLengths = SummaryFeatures.GetSegLengths(ObservedCNPs);
    }

    public void Optimize(SimParams simParams, Random rnd, int repeats)
    {
        GenerateSimulatedCNPs(simParams, rnd, repeats);
        var segDist = GetSegLengthDistance();
        var cpDist = GetChangepointDistance();
        var bpDist = GetBreakpointDistance();
        Console.WriteLine($"Seg Length WD: {segDist}; Changepoint WD: {cpDist}; BP per chr WD: {bpDist}");
        return;
    }
    public void GenerateSimulatedCNPs(SimParams simParams, Random rnd, int repeats)
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
        SimulatedCNPs = GetCNPs(GenRef, samples);
    }
    public double GetSegLengthDistance()
    {
        var dataSegList = SummaryFeatures.GetSegLengths(ObservedCNPs);
        var simSegList = SummaryFeatures.GetSegLengths(SimulatedCNPs);
        var histMax = GenRef.ChrLengths["chr1"];
        var histMin = 0;
        var histBins = 1000;
        return CalculateDistance(dataSegList, simSegList, histBins, histMin, histMax);
    }

    public double CalculateDistance(List<double> data, List<double> sim, int bins, int min, int max)
    {
        var dataHist = new Histogram(data, bins, min, max);
        var simHist  = new Histogram(sim, bins, min, max);
        return StatisticMeasures.WassersteinDistance(dataHist, simHist);
    }

    public double GetChangepointDistance()
    {
        var dataChangeInfo = SummaryFeatures.GetChangepointInfo(ObservedCNPs);
        var simChangeInfo  = SummaryFeatures.GetChangepointInfo(SimulatedCNPs);
        var histMax = Math.Max(dataChangeInfo.max, simChangeInfo.max);
        var histMin = 0;
        var histBins = 50;
        return CalculateDistance(dataChangeInfo.values, simChangeInfo.values, histBins, histMin, histMax);

    }

    public double GetBreakpointDistance()
    {
        var dataBPInfo = SummaryFeatures.GetBreakpointsPerChromosome(ObservedCNPs);
        var simBPInfo  = SummaryFeatures.GetBreakpointsPerChromosome(SimulatedCNPs);
        var histMax = Math.Max(dataBPInfo.max, simBPInfo.max);
        var histMin = 0;
        var histBins = 50;
        return CalculateDistance(dataBPInfo.values, simBPInfo.values, histBins, histMin, histMax);
    }


    public void SetSimulatedDistribution(List<Sample> simulatedData)
        => SimulatedCNPs = GetCNPs(GenRef, simulatedData);

    public double GetTotalDistance()
    {
        //var segDist = StatisticMeasures<long>.WassersteinDistance(ObservedSegLengths, SimulatedSegLengths);
        //var cpDist = StatisticMeasures<int>.WassersteinDistance(ObservedChangepoints, SimulatedChangepoints);
        //var bpDist = StatisticMeasures<long>.WassersteinDistance(ObservedBreakpoints, SimulatedBreakpoints);
        return 0;//segDist*segDist + cpDist*cpDist + bpDist*bpDist;
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