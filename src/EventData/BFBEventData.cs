using SimChA.Data;
using SimChA.Simulation;

namespace SimChA.EventData;

public record BFBEventData : ContigEventData
{
    public bool Direction { get; }
    public long Start { get; }
    public long FinalBreak { get; }

    public BFBEventData(
        Random rnd,
        CNEventPars cnEventPars,
        int contigId,
        long contigLen,
        TerminalArm arm)
        : base(cnEventPars, contigId, contigLen)
    {
        if (arm.UsableLength <= 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(arm), "BFB requires a non-empty terminal arm.");
        }

        // The first break removes a terminal segment but cannot enter the centromere. Frac keeps
        // its historical contig-wide scale; an overlong draw collapses to a whole-arm break.
        long drawnLength = Sampling.GetExpSeg(rnd, contigLen, cnEventPars.Frac);
        long removedLength = Math.Min(drawnLength, arm.UsableLength);
        Direction = arm.Direction;
        Start = Direction ? removedLength : contigLen - removedLength;

        // After sister fusion the two centromeres flank a symmetric copy of the retained terminal
        // arm. A uniformly placed cut in that interval finishes the BFB and leaves one centromere
        // in the retained daughter.
        long gapStart = contigLen - arm.UsableLength;
        long gapEnd = contigLen + arm.UsableLength - 2 * removedLength;
        FinalBreak = gapStart == gapEnd
            ? gapStart
            : rnd.NextInt64(gapStart, gapEnd + 1);
    }

    public override void ApplyEvent(Karyotype kar)
        => kar.ApplyBFB(ContigId, Start, Direction, FinalBreak);

    public override string EventDesc()
        => base.EventDesc() +
           $"start:{(Direction ? 0 : Start)};end:{(Direction ? Start : Length)};" +
           $"final_break:{FinalBreak};";
}
