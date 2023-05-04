using SimChA.DataTypes;
using SimChA.Misc;
using SimChA.Simulation;

namespace SimChA.EventData;

public record TailEventData : ContigEventData
{
    public readonly long DelFraction= -1;
    public readonly bool Direction;
    
    // Constructor used for Tail Events
    public TailEventData(Random rnd, Karyotype kar, CNEventP eventP, int contigId) : base(eventP, contigId)
    {
        long tailSize = eventP.Get("Size", 1_000_000);
        DelFraction = Sampling.GetExpSeg(rnd, kar.ContigLen(contigId), tailSize);
        Direction = rnd.CoinFlip();
    }

    public override string ToString()
    {
        return EventType switch
        {
            CNEventType.TailDeletion => $"{EventType}\t{ContigId}\t{DelFraction}\t{Direction}",
            CNEventType.BreakageFusionBridge => $"{EventType}\t{ContigId}\t{DelFraction}\t{Direction}",
            _ => throw new ArgumentOutOfRangeException(nameof(EventType), EventType, null)
        };
    }
    
    public override string ApplyEvent(Karyotype kar)
        => kar.ApplyEvent(this);
}
