// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using Extreme.Mathematics;
using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.EventData;

public record PyrgoEventData : BaseEventData
{
    public int ContigId { get; }
    public List<(long, long)> FragmentsList { get; }
    
    private static List<(long, long)> GetSubsegments(Random rnd, long start, long fragmentLen, double mean)
    {
        var meanSize = (long) (fragmentLen * mean);
        int fracCount = Sampling.GetFragCount(rnd, mean);
        var frags = new List<(long, long)>();
        for (int i = 0; i < fracCount; i++)
        {
            long fracLen = Sampling.GetExpSeg(rnd, fragmentLen, meanSize);
            long fracStart = Sampling.GetInternalPos(rnd, fragmentLen - fracLen);
            frags.Add((start + fracStart, fracLen));
        }
        return frags;
    }
    
    public PyrgoEventData(Random rnd, Karyotype kar, CNEventP eventP, int contigId) : base(eventP)
    {
        ContigId = contigId;
        long contigLen = kar.ContigLen(contigId);
        long pyrgoLen = eventP.Get("Size", 1_000_000L);
        double pyrgoMean = eventP.Get("Mean", 0.1);
        long pyrgoFrag = Sampling.GetExpSeg(rnd, contigLen, pyrgoLen);
        long pyrgoStart = Sampling.GetInternalPos(rnd, contigLen - pyrgoFrag);
        FragmentsList = GetSubsegments(rnd, pyrgoStart, pyrgoFrag, pyrgoMean);
    }
    
    public override string ToString() => $"{EventType}\t{ContigId}";
    
    public override string ApplyEvent(Karyotype kar)
        => kar.ApplyEvent(this);
}