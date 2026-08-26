using SimChA.Computation;
using SimChA.Data;
using SimChA.Simulation;

namespace SimChA.EventData;

public record InternalEventData : ContigEventData
{
    public bool Direction { get; }
    public long Start { get; }
    public long End { get; }

    // Constructor used for internal events
    public InternalEventData(Random rnd, CNEventPars CNEventPars, int contigId, long contigLen)
        : this(rnd, CNEventPars, contigId, contigLen, true, contigLen)
    { }

    // Constructor used after a terminal chromosome arm has been selected. Event length is drawn
    // first relative to that arm, then the start is uniform over every placement that stays within
    // the arm and therefore cannot cross its inward centromere boundary.
    public InternalEventData(
        Random rnd,
        CNEventPars CNEventPars,
        int contigId,
        long contigLen,
        bool direction,
        long armLength) : base(CNEventPars, contigId, contigLen)
    {
        Direction = direction;
        long segLen = Sampling.GetParetoSeg(rnd, armLength, CNEventPars.Frac);
        long armStart = direction ? 0 : contigLen - armLength;
        long lastStart = armStart + armLength - segLen;
        Start = lastStart == armStart
            ? armStart
            : rnd.NextInt64(armStart, lastStart + 1);
        End = Start + segLen;
    }

    // Constructor for Centromere-bound events
    public InternalEventData(Random rnd, CNEventPars CNEventPars, int contigId, IEnumerable<(long start, long end)> centromeres, long contigLen) 
        : base(CNEventPars, contigId, contigLen)
    {
        (long start, long end) = centromeres.Shuffle(rnd).First();
        long pos = rnd.NextInt64(start, end);
        long segLen = Sampling.GetExpSeg(rnd, contigLen, CNEventPars.Frac);
        if (rnd.CoinFlip())
        {
            Start = pos;
            End = Math.Min(pos + segLen, contigLen);
        }
        else
        {
            Start = Math.Max(0, pos - segLen);
            End = pos;
        }
    }

    public override void ApplyEvent(Karyotype kar)
    {
        switch (EventType)
        {
            case CNEventType.InternalDuplication:
            case CNEventType.CentromereBoundDuplication:
                kar.ApplyInternalDuplication(ContigId, Start, End);
                break;
            case CNEventType.InternalDeletion:
            case CNEventType.CentromereBoundDeletion:
                kar.ApplyInternalDeletion(ContigId, Start, End);
                break;
            case CNEventType.InternalInversion:
                kar.ApplyInternalInversion(ContigId, Start, End);
                break;
            case CNEventType.InvertedDuplication:
                kar.ApplyInvertedDuplication(ContigId, Start, End);
                break;
            default:
                throw new Exception($"Invalid event type {EventType} for InternalEventData");
        }
    }

    public override string EventDesc()
        => base.EventDesc() + $"start:{Start};end:{End}";
}
