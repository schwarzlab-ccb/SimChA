using SimChA.Computation;
using SimChA.Simulation;

namespace SimChA.EventData;

public record TailEventData : ContigEventData
{
    public long DelFraction { get; }
    public bool Direction { get; }
    
    // Constructor used for Tail CNEventPars
    public TailEventData(Random rnd, CNEventPars CNEventPars, int contigId, long contigLen) : base(CNEventPars, contigId)
    {
        long tailSize = CNEventPars.GetInt("Size");
        DelFraction = Sampling.GetExpSeg(rnd, contigLen, tailSize);
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
