// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.Simulation;

namespace SimChA.EventData;

public record RigmaEventData : ContigEventData
{
    public long Start { get; }
    public List<long> StopsList { get; }
    
    public RigmaEventData(Random rnd, CNEventPars cnEventPars, int contigId, long contigLen) : base(cnEventPars, contigId)
    {
        ContigId = contigId;
        long rigmaLen = cnEventPars.Get("Size", 1_000_000L);
        double rigmaMean = cnEventPars.Get("Mean", 0.1);
        Start= Sampling.GetInternalPos(rnd, contigLen - rigmaLen);
        int rigmaCount = Sampling.GetFragCount(rnd, rigmaMean);
        StopsList = Enumerable.Range(0, rigmaCount).Select(_ => Sampling.GetExpSeg(rnd, contigLen, rigmaMean)).ToList();
    }
    
    public override void ApplyEvent(Karyotype kar)
        => kar.ApplyRigma(ContigId, Start, StopsList);

    public override string ToString()
        => $"contig:{ContigId};start{Start};stops:{string.Join(",", StopsList)}";
}
