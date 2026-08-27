using SimChA.Computation;
using SimChA.Data;
using SimChA.Simulation;

namespace SimChA.EventData;

public record InternalEventData : ContigEventData
{
    public long Start { get; }
    public long End { get; }

    // Constructor used after a terminal chromosome arm has been selected. Interior event lengths
    // follow a Pareto (shape 0.5) fit to the empirical proportion distribution, drawn against the
    // fitted arm scale and then bounded by the arm proper. The start is uniform over every
    // placement that keeps the event inside that arm, so it never reaches the centromere.
    public InternalEventData(
        Random rnd,
        CNEventPars CNEventPars,
        int contigId,
        long contigLen,
        TerminalArm arm) : base(CNEventPars, contigId, contigLen)
    {
        long segLen = Math.Min(
            Sampling.GetParetoSeg(rnd, arm.ArmLength, CNEventPars.Frac), arm.UsableLength);
        long armStart = arm.Direction ? 0 : contigLen - arm.UsableLength;
        Start = rnd.NextInt64(armStart, armStart + arm.UsableLength - segLen + 1);
        End = Start + segLen;
    }

    // Constructor for Centromere-bound events, whose lengths stay exponential
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
