using SimChA.Computation;
using SimChA.Data;

namespace SimChA.EventData;

public record PairEventData : BaseEventData
{
    public readonly int ContigIdA;
    public readonly int ContigIdB;
    public readonly long PosA;
    public readonly long PosB;
    public readonly bool Inverted;

    private static long SampleArmBreakpoint(Random rnd, long contigLength, TerminalArm arm)
    {
        if (arm.UsableLength <= 0 || arm.UsableLength > contigLength)
        {
            throw new ArgumentOutOfRangeException(
                nameof(arm), "Translocation requires a non-empty terminal arm.");
        }
        return arm.Direction
            ? rnd.NextInt64(0, arm.UsableLength)
            : rnd.NextInt64(contigLength - arm.UsableLength, contigLength);
    }

    // Breakpoints are sampled outside the centromeres. If the selected arms face opposite
    // directions, B is reversed and its breakpoint is converted into reversed-contig coordinates;
    // the exchanged tails then leave exactly one centromere on each derivative.
    public PairEventData(
        Random rnd,
        CNEventPars cnEventPars,
        int contigA,
        long lenA,
        TerminalArm armA,
        int contigB,
        long lenB,
        TerminalArm armB)
        : base(cnEventPars)
    {
        ContigIdA = contigA;
        ContigIdB = contigB;
        PosA = SampleArmBreakpoint(rnd, lenA, armA);
        long originalPosB = SampleArmBreakpoint(rnd, lenB, armB);
        Inverted = armA.Direction != armB.Direction;
        PosB = Inverted ? lenB - originalPosB : originalPosB;
    }
    
    public override void ApplyEvent(Karyotype kar)
        => kar.ApplyTranslocation(ContigIdA, ContigIdB, PosA, PosB, Inverted);
    
    public override string EventDesc()
        => $"contigA:{ContigIdA};contigB:{ContigIdB};posA:{PosA};posB:{PosB};invertedB:{Inverted}";
}