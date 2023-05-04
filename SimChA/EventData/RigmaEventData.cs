// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.EventData;

public record RigmaEventData : BaseEventData
{
    public int ContigId { get; }
    public long Start { get; }
    public List<long> StopsList { get; }
    
    public RigmaEventData(Random rnd, Karyotype kar, CNEventP eventP, int contigId) : base(eventP)
    {
        ContigId = contigId;
        long contigLen = kar.ContigLen(contigId);
        long rigmaLen = eventP.Get("Size", 1_000_000L);
        double rigmaMean = eventP.Get("Mean", 0.1);
        Start= Sampling.GetInternalPos(rnd, contigLen - rigmaLen);
        int rigmaCount = Sampling.GetFragCount(rnd, rigmaMean);
        StopsList = Enumerable.Range(0, rigmaCount).Select(i => Sampling.GetExpSeg(rnd, contigLen, rigmaMean)).ToList();
    }
    
    public override void ApplyEvent(Karyotype kar)
        => kar.ApplyRigma(ContigId, Start, StopsList);

    public override string ToString()
        => $"contig:{ContigId};start{Start};stops:{string.Join(",", StopsList)}";
}
