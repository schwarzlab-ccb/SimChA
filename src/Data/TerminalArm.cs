namespace SimChA.Data;

// One physical contig end together with the arm reaching inward from it to the nearest centromere.
// `ArmLength` runs to the middle of that centromere: it is the length scale the empirical event
// proportions are fitted against, because arms are cut at the centromere midpoint upstream.
// `UsableLength` stops at the near edge of the centromere and bounds what an event may actually
// cover, so a drawn length longer than the arm proper collapses onto a whole-arm event instead of
// entering the centromere. `Direction` is true for the front/5' end and false for the back/3' end.
public readonly record struct TerminalArm(bool Direction, long ArmLength, long UsableLength);
