// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using Extreme.Mathematics;
using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.EventData;

public record PyrgoEventData : BaseEventData
{
    public int ContigId { get; }
    public List<(long, long)> FragmentsList { get; }
    
    public PyrgoEventData(Random rnd, Karyotype kar, CNEventP eventP, int contigId) : base(eventP)
    {
        ContigId = contigId;
        long contigLen = kar.ContigLen(contigId);
        long pyrgoLen = eventP.Get("Size", 1_000_000L);
        double pyrgoMean = eventP.Get("Mean", 0.1);
        long pyrgoFrag = Sampling.GetExpSeg(rnd, contigLen, pyrgoLen);
        long pyrgoStart = Sampling.GetInternalPos(rnd, contigLen - pyrgoFrag);
        
        var meanSize = (long)(pyrgoFrag * pyrgoMean);
        int fracCount = Sampling.GetFragCount(rnd, pyrgoMean);
        FragmentsList = new List<(long, long)>();
        for (int i = 0; i < fracCount; i++)
        {
            long fracLen = Sampling.GetExpSeg(rnd, pyrgoFrag, meanSize);
            long fracStart = Sampling.GetInternalPos(rnd, pyrgoFrag - fracLen);
            FragmentsList.Add((pyrgoStart + fracStart, fracLen));
        }
    }

    public override void ApplyEvent(Karyotype kar)
        => kar.ApplyPyrgo(ContigId, FragmentsList);

    public override string ToString()
        => $"contig:{ContigId};frags:{string.Join(",", FragmentsList)}";
}