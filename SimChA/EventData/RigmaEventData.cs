using Extreme.Statistics.Distributions;
using SimChA.Data;
using SimChA.Simulation;

namespace SimChA.EventData;

public record RigmaEventData : ContigEventData
{
    public long Start { get; }
    public List<long> StopsList { get; }
    
    public RigmaEventData(Random rnd, CNEventPars cnEventPars, int contigId, long contigLen) : base(cnEventPars, contigId)
    {
        ContigId = contigId;
        int fracCount = GeometricDistribution.Sample(rnd, 1 / cnEventPars.Frag) + 1;
        Start = Sampling.GetPos(rnd, contigLen - cnEventPars.Size);
        StopsList = Enumerable.Range(0, fracCount).Select(_ => Sampling.GetExpSeg(rnd, contigLen, cnEventPars.Frag / cnEventPars.Size)).ToList();
    }
    
    public override void ApplyEvent(Karyotype kar)
        => kar.ApplyRigma(ContigId, Start, StopsList);

    public override string ToString()
        => $"contig:{ContigId};start{Start};stops:{string.Join(",", StopsList)}";
}
