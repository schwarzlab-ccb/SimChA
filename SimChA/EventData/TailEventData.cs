using System.Runtime.InteropServices;
using SimChA.Computation;
using SimChA.Simulation;

namespace SimChA.EventData;

public record TailEventData : ContigEventData
{
    public long Length { get; }
    public bool Direction { get; }
    
    // Constructor used for Tail CNEventPars
    public TailEventData(Random rnd, CNEventPars CNEventPars, int contigId, long contigLen) : base(CNEventPars, contigId)
    {
        Length = Sampling.GetPos(rnd, contigLen);
        Direction = rnd.CoinFlip();
    }

    public override void ApplyEvent(Karyotype kar)
    {
        if (EventType == CNEventType.TailDeletion)
        {
            kar.ApplyTailDeletion(ContigId, Length, Direction);
        }
        else if (EventType == CNEventType.TailDuplication)
        {
            kar.ApplyTailDuplication(ContigId, Length, Direction);
        }
        else if (EventType == CNEventType.BreakageFusionBridge)
        {
            kar.ApplyBFB(ContigId, Length, Direction);
        }
    }
    
    public override string ToString()
        => $"contig:{ContigId};length:{Length};dir:{Direction}";
}
