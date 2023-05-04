using SimChA.DataTypes;
using SimChA.Misc;
using SimChA.Simulation;

namespace SimChA.EventData;

public record TailEventData : ContigEventData
{
    public long DelFraction { get; }
    public bool Direction { get; }
    
    // Constructor used for Tail Events
    public TailEventData(Random rnd, Karyotype kar, CNEventP eventP, int contigId) : base(eventP, contigId)
    {
        long tailSize = eventP.Get("Size", 1_000_000);
        DelFraction = Sampling.GetExpSeg(rnd, kar.ContigLen(contigId), tailSize);
        Direction = rnd.CoinFlip();
    }

    public override void ApplyEvent(Karyotype kar)
    {
        if (EventType == CNEventType.TailDeletion)
        {
            kar.ApplyTailDeletion(ContigId, DelFraction, Direction);
        }
        else if (EventType == CNEventType.BreakageFusionBridge)
        {
            kar.ApplyBFB(ContigId, DelFraction, Direction);
        }
    }
    
    public override string ToString()
        => $"contig:{ContigId};delFraction:{DelFraction};dir:{Direction}";
}
