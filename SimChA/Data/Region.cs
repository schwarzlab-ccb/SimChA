namespace SimChA.Data;

// A region is zero indexed start-inclusive, end-exclusive, e.g. [0, 1) is a region of length 1 containing the first base.
// TODO: The direction should not be part of the region, should be part of the chromosome
public record Region(long Start, long End, string ChrNo, bool Hap1, bool Forward = true, 
                        List<SNV>? SNVs = null) : GenRange(Start, End, ChrNo)
{
    private static string HapToString(bool parent) 
        => parent ? "H1" : "H2";

    private string HapString => HapToString(Hap1);

    public static string DirToStr(bool dir) => dir ? ">" : "<";

    private string DirString => DirToStr(Forward);
    
    public override string ToString() => $"{HapString}{DirString}{ChrNo}[{Start}:{End})";

    public int NumSNVsBetween(long start, long end)
        => SNVs != null && SNVs.Count > 0 ? SNVs.Count(s => s.Location >= start && s.Location <= end) : 0;
}

public record PArm(long Start, long End, string ChrNo, bool Hap1, bool Forward = true, List<SNV>? SNVs = null) 
    : Region(Start, End, ChrNo, Hap1, Forward, SNVs)
{
    public override string ToString() => base.ToString();
}

public record QArm(long Start, long End, string ChrNo, bool Hap1, bool Forward = true, List<SNV>? SNVs = null) 
    : Region(Start, End, ChrNo, Hap1, Forward, SNVs)
{
    public override string ToString() => base.ToString();
}

public record Centromere(long Start, long End, string ChrNo, bool Hap1, bool Forward = true, List<SNV>? SNVs = null) 
    : Region(Start, End, ChrNo, Hap1, Forward, SNVs)
{ 
    public override string ToString() => base.ToString();
}