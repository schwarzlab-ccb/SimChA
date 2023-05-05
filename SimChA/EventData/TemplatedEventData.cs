// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using Extreme.Statistics.Distributions;
using SimChA.DataTypes;
using SimChA.Misc;
using SimChA.Simulation;

namespace SimChA.EventData;

public record TemplatedEventData : BaseEventData
{
    public List<(int id, long start, long len, bool dir)> Frags { get; } = new();
    
    public TemplatedEventData(Random rnd, Karyotype kar, CNEventP eventP, IEnumerable<int> seq) : base(eventP)
    {
        var size = eventP.Get("Size", 1_000_000L);
        var fragMean = eventP.Get("Frag", 10.0);
        var contigCount = GeometricDistribution.Sample(rnd, 1 / fragMean) 
                          + (eventP.Type != CNEventType.TIBridge ? 1 : 2);
        
        var contigIds = seq.Take(contigCount).ToList();
        for (var i = 0; i < contigIds.Count; i++)
        {
            var id = contigIds[i];
            var contigLen = kar.ContigLen(id); 
            // First segment of a bridge, or first and last on a chain do not have a length
            bool skipLen = i == 0 && eventP.Type != CNEventType.TICycle ||
                           i == contigCount - 1 && eventP.Type == CNEventType.TIChain;
            var fragLen = skipLen ? 0L : Sampling.GetExpSeg(rnd, contigLen, size);
            var fragStart = Sampling.GetInternalPos(rnd, contigLen - fragLen);
            var dir = i == 0 || rnd.CoinFlip();
            Frags.Add((id, fragStart, fragLen, dir));
        }
    }

    public override void ApplyEvent(Karyotype kar)
    {
        if (EventType == CNEventType.TIBridge)
        {
            kar.ApplyTIBridge(Frags);
        }
        else if (EventType == CNEventType.TIChain)
        {
            kar.ApplyTIChain(Frags);
        }
        else if (EventType == CNEventType.TICycle)
        {
            kar.ApplyTICycle(Frags);
        }
    }

    public override string ToString()
        => string.Join(",", Frags.Select(x => $"({x.id},{x.start},{x.len},{Region.DirToStr(x.dir)})"));
}