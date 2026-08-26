using SimChA.Computation;
using SimChA.Data;
using SimChA.Simulation;

namespace SimChA.EventData;

public record TailEventData : ContigEventData
{
    public bool Direction { get; }
    public long Start { get; }
    
    // Constructor used for Tail CNEventPars
    public TailEventData(Random rnd, CNEventPars CNEventPars, int contigId, long contigLen) 
        : this(CNEventPars, contigId, ComputeTailParams(rnd, CNEventPars, contigLen), contigLen)
    { }

    // Constructor used for events constrained to a preselected physical contig end.
    public TailEventData(
        Random rnd,
        CNEventPars CNEventPars,
        int contigId,
        long contigLen,
        bool direction,
        long? maxLength = null)
        : base(CNEventPars, contigId, contigLen)
    {
        long segLen = Sampling.GetExpSeg(rnd, contigLen, CNEventPars.Frac);
        if (maxLength is { } limit)
        {
            segLen = Math.Min(segLen, limit);
        }
        Direction = direction;
        Start = Direction ? segLen : contigLen - segLen;
    }
    
    public TailEventData(Random rnd, CNEventPars CNEventPars, int contigId, IEnumerable<(long start, long end)> cents, long contigLen)
        : this(CNEventPars, contigId, ComputeArmParams(rnd, cents, contigLen), contigLen)
    { }
    
    private TailEventData(CNEventPars CNEventPars, int contigId, (long segLen, bool direction) segment, long contigLen)
        : base(CNEventPars, contigId, contigLen)
    {
        Direction = segment.direction;
        Start = Direction ? segment.segLen : contigLen - segment.segLen;
    }

    private static (long segLen, bool direction) ComputeTailParams(
        Random rnd,
        CNEventPars CNEventPars,
        long contigLen)
    {
        // Preserve the historical seeded-simulation order: length is sampled before direction.
        long segLen = Sampling.GetExpSeg(rnd, contigLen, CNEventPars.Frac);
        bool direction = rnd.CoinFlip();
        return (segLen, direction);
    }
    
    private static (long segLen, bool direction) ComputeArmParams(
        Random rnd, IEnumerable<(long start, long end)> cents, long contigLen)
    {
        var cent = cents.Shuffle(rnd).First();
        bool direction = rnd.CoinFlip();
        long length = direction ? cent.start : contigLen - cent.end;
        long breakpoint = Math.Max(1, rnd.NextInt64(cent.end - cent.start)); // Somewhere within the centromere, uniform 
        length += breakpoint;
        return (length, direction);
    }

    public override void ApplyEvent(Karyotype kar)
    {
        switch (EventType)
        {
            case CNEventType.TailDeletion:
            case CNEventType.TelomereDeletion:
            case CNEventType.ArmDeletion:
                kar.ApplyTailDeletion(ContigId, Start, Direction);
                break;
            case CNEventType.TailDuplication:
            case CNEventType.TelomereDuplication:
                kar.ApplyTailDuplication(ContigId, Start, Direction);
                break;
            case CNEventType.ArmDuplication:
                kar.ApplyDetachedTailDuplication(ContigId, Start, Direction);
                break;
            case CNEventType.BreakageFusionBridge:
                kar.ApplyBFB(ContigId, Start, Direction);
                break;
            default:
                throw new Exception($"Invalid event type {EventType} for TailEventData");
        }
    }
    
    public override string EventDesc()
        => base.EventDesc() + $"start:{(Direction ? 0 : Start)};end:{(Direction ? Start : Length)};";
}
