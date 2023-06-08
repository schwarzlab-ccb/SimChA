// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using Extreme.Statistics.Distributions;
using SimChA.Simulation;

namespace SimChA.EventData;

public record RigmaEventData : ContigEventData
{
    public long Start { get; }
    public List<long> StopsList { get; }
    
    public RigmaEventData(Random rnd, CNEventPars cnEventPars, int contigId, long contigLen) : base(cnEventPars, contigId)
    {
        ContigId = contigId;
        long rigmaLen = cnEventPars.GetLong("Size");
        double fragMean = cnEventPars.GetDouble("Frag");
        int fracCount = GeometricDistribution.Sample(rnd, 1 / fragMean) + 1;
        Start = Sampling.GetInternalPos(rnd, contigLen - rigmaLen);
        StopsList = Enumerable.Range(0, fracCount).Select(_ => Sampling.GetExpSeg(rnd, contigLen, fragMean / rigmaLen)).ToList();
    }
    
    public override void ApplyEvent(Karyotype kar)
        => kar.ApplyRigma(ContigId, Start, StopsList);

    public override string ToString()
        => $"contig:{ContigId};start{Start};stops:{string.Join(",", StopsList)}";
}
