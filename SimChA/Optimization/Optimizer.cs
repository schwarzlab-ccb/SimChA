using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.Simulation;
namespace SimChA.Optimization;

public class Optimizer
{   
    public Dictionary<string, List<CopyNumber>> ObservedCNPs { get; }
    public List<long> ObservedSegLengths { get; }
    public List<int> ObservedChangepoints { get; }
    public List<long> ObservedBreakpoints { get; }

    public Dictionary<string, List<CopyNumber>> SimulatedCNPs { get; set;}
    public List<long> SimulatedSegLengths { get; set;}
    public List<int> SimulatedChangepoints { get; set;}
    public List<long> SimulatedBreakpoints { get; set;}
    GenRef GenRef { get; }
    public Optimizer(GenRef genRef, List<Sample> observedData)
    {
        GenRef = genRef;
        ObservedCNPs = GetCNPs(GenRef, observedData);
        ObservedSegLengths = SummaryFeatures.GetSegLengths(ObservedCNPs);
        ObservedChangepoints = SummaryFeatures.GetChangepoints(ObservedCNPs);
        ObservedBreakpoints = SummaryFeatures.GetBreakpointsPerChromosome(ObservedCNPs);
    }

    public void SetSimulatedDistribution(List<Sample> simulatedData)
    {
        SimulatedCNPs = GetCNPs(GenRef, simulatedData);
        SimulatedSegLengths = SummaryFeatures.GetSegLengths(SimulatedCNPs);
        SimulatedChangepoints = SummaryFeatures.GetChangepoints(SimulatedCNPs);
        SimulatedBreakpoints = SummaryFeatures.GetBreakpointsPerChromosome(SimulatedCNPs);
    }

    public double GetTotalDistance()
    {
        var segDist = StatisticMeasures<long>.WassersteinDistance(ObservedSegLengths, SimulatedSegLengths);
        var cpDist = StatisticMeasures<int>.WassersteinDistance(ObservedChangepoints, SimulatedChangepoints);
        var bpDist = StatisticMeasures<long>.WassersteinDistance(ObservedBreakpoints, SimulatedBreakpoints);
        return segDist + cpDist + bpDist;
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