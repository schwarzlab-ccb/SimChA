// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using Extreme.Statistics.Distributions;
using SimChA.DataTypes;
using SimChA.Misc;
using SimChA.Simulation;

namespace SimChA.EventData;

public record TemplatedEventData : BaseEventData
{
    public List<(int, long, long, bool)> Frags { get; }
    
    public TemplatedEventData(Random rnd, Karyotype kar, CNEventP eventP, IEnumerator<int> IDsEnumerator) : base(eventP)
    {
        var size = eventP.Get("Size", 1_000_000L);
        var fragMean = eventP.Get("Frag", 10.0);
        var fragCount = GeometricDistribution.Sample(rnd, 1 / fragMean) 
                        + (eventP.Type != CNEventType.TIBridge ? 1 : 2);
        var fragments = new List<(int id, long start, long len, bool dir)>();
        for (var i = 0; i < fragCount; i++, IDsEnumerator.MoveNext())
        {
            var id = IDsEnumerator.Current;
            var contigLen = kar.ContigLen(id); 
            // First segment of a bridge, or first and last on a chain do not have a length
            bool skipLen = i == 0 && eventP.Type != CNEventType.TICycle ||
                           i == fragCount - 1 && eventP.Type == CNEventType.TIChain;
            var fragLen = skipLen ? 0L : Sampling.GetExpSeg(rnd, contigLen, size);
            var fragStart = Sampling.GetInternalPos(rnd, contigLen - fragLen);
            var dir = i == 0 || rnd.CoinFlip();
            fragments.Add((id, fragStart, fragLen, dir));
        }
        Frags = fragments;
    }
    
    public override string ToString() => $"{EventType}\t{string.Join(",", Frags)}";
    
    public override string ApplyEvent(Karyotype kar)
        => kar.ApplyEvent(this);
}