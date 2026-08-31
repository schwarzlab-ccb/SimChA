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

    // Constructor used after a terminal chromosome arm has been selected. Telomere-bound lengths
    // follow the fitted Beta, the remaining terminal events stay exponential; both are drawn
    // against the fitted arm scale and then bounded by the arm proper, so a draw longer than the
    // arm becomes a whole-arm event rather than one reaching into the centromere.
    public TailEventData(
        Random rnd,
        CNEventPars CNEventPars,
        int contigId,
        long contigLen,
        TerminalArm arm)
        : base(CNEventPars, contigId, contigLen)
    {
        long drawnLen = CNEventPars.Type is CNEventType.TelomereDeletion or CNEventType.TelomereDuplication
            ? Sampling.GetBetaSeg(rnd, arm.ArmLength, CNEventPars.Frac)
            : Sampling.GetExpSeg(rnd, arm.ArmLength, CNEventPars.Frac);
        long segLen = Math.Min(drawnLen, arm.UsableLength);
        Direction = arm.Direction;
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
            default:
                throw new Exception($"Invalid event type {EventType} for TailEventData");
        }
    }
    
    public override string EventDesc()
        => base.EventDesc() + $"start:{(Direction ? 0 : Start)};end:{(Direction ? Start : Length)};";
}
