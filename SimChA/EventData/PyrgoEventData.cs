// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using Extreme.Statistics.Distributions;
using SimChA.Simulation;

namespace SimChA.EventData;

public record PyrgoEventData : ContigEventData
{
    public List<(long start, long length)> FragmentsList { get; }
    
    public PyrgoEventData(Random rnd, CNEventPars cnEventPars, int contigId, long contigLen) : base(cnEventPars, contigId)
    {
        ContigId = contigId;
        long pyrgoFrag = Sampling.GetExpSeg(rnd, contigLen, cnEventPars.Size);
        long pyrgoStart = Sampling.GetInternalPos(rnd, contigLen - pyrgoFrag);
        
        int fracCount = GeometricDistribution.Sample(rnd, 1.0 / cnEventPars.Frag) + 1;
        
        FragmentsList = new List<(long, long)>();
        for (int i = 0; i < fracCount; i++)
        {
            long fracLen = Sampling.GetExpSeg(rnd, pyrgoFrag, cnEventPars.Size / cnEventPars.Frag);
            long fracStart = Sampling.GetInternalPos(rnd, pyrgoFrag - fracLen);
            FragmentsList.Add((pyrgoStart + fracStart, fracLen));
        }
    }

    public override void ApplyEvent(Karyotype kar)
        => kar.ApplyPyrgo(ContigId, FragmentsList);

    public override string ToString()
        => $"contig:{ContigId};frags:{string.Join(",", FragmentsList)}";
}