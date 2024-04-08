using System.Runtime.InteropServices;
using SimChA.Computation;
using SimChA.Simulation;

namespace SimChA.EventData;

public record ArmEventData : ContigEventData
{
    public long Fraction { get; }
    public bool Direction { get; }
    
    // Constructor used for Tail CNEventPars
    public ArmEventData(Random rnd, CNEventPars CNEventPars, int contigId, long contigLen) : base(CNEventPars, contigId)
    {
        Fraction = Sampling.GetPos(rnd, contigLen);
        Direction = rnd.CoinFlip();
    }

    public override void ApplyEvent(Karyotype kar)
    {
        if (EventType == CNEventType.TailDeletion)
        {
            kar.ApplyTailDeletion(ContigId, Fraction, Direction);
        }
        else if (EventType == CNEventType.TailDuplication)
        {
            kar.ApplyTailDuplication(ContigId, Fraction, Direction);
        }
        else if (EventType == CNEventType.BreakageFusionBridge)
        {
            kar.ApplyBFB(ContigId, Fraction, Direction);
        }
    }
    
    public override string ToString()
        => $"contig:{ContigId};fraction:{Fraction};dir:{Direction}";
}
