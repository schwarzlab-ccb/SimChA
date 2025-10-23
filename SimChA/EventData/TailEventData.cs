using SimChA.Computation;
using SimChA.Data;
using SimChA.Simulation;

namespace SimChA.EventData;

public record TailEventData : ContigEventData
{
    public bool Direction { get; }
    private long Start { get; }
    
    // Constructor used for Tail CNEventPars
    public TailEventData(Random rnd, CNEventPars CNEventPars, int contigId, long contigLen) 
        : base(CNEventPars, contigId, contigLen)
    {
        Length = Sampling.GetExpSeg(rnd, contigLen, CNEventPars.Frac);
        Direction = rnd.CoinFlip();
        Start = Direction ? 0 : contigLen - Length;
    }
    
    public TailEventData(Random rnd, CNEventPars CNEventPars, int contigId, IEnumerable<(long start, long end)> cents, long contigLen)
        : base(CNEventPars, contigId, -1)
    {
        var cent =  cents.Shuffle(rnd).First();
        Direction = rnd.CoinFlip();
        Length = Direction ? cent.start : contigLen - cent.end;
        long breakpoint = Math.Max(1, rnd.NextInt64(cent.end - cent.start)); // Somewhere within the centromere, uniform 
        Length += breakpoint;
        Start = Direction ? 0 : contigLen - Length;
    }

    public override void ApplyEvent(Karyotype kar)
    {
        switch (EventType)
        {
            case CNEventType.TailDeletion:
            case CNEventType.ArmDeletion:
                kar.ApplyTailDeletion(ContigId, Length, Direction);
                break;
            case CNEventType.TailDuplication:
            case CNEventType.ArmDuplication:
                kar.ApplyTailDuplication(ContigId, Length, Direction);
                break;
            case CNEventType.BreakageFusionBridge:
                kar.ApplyBFB(ContigId, Length, Direction);
                break;
            default:
                throw new Exception($"Invalid event type {EventType} for TailEventData");
        }
    }
    
    public override string EventDesc()
        => base.EventDesc() +  $"start:{Start};end:{Start + Length};";
}
