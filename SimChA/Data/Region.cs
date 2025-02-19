using System.Text;

namespace SimChA.Data;

// A region is zero indexed start-inclusive, end-exclusive, e.g. [0, 1) is a region of length 1 containing the first base.
public record Region(long Start, long End, string ChrNo, bool Hap1, List<SNV>? SNVs = null) : GenRange(Start, End, ChrNo)
{
    private static string HapToString(bool parent) 
        => parent ? "H1" : "H2";

    private string HapString 
        => HapToString(Hap1);
    
    public override string ToString() 
        => $"{HapString}{base.ToString()}";

    public int NumSNVsBetween(long start, long end)
        => SNVs is { Count: > 0 } ? SNVs.Count(s => s.Pos >= start && s.Pos <= end) : 0;

    public string GetSeq(GenRef genRef)
    {
        char[] regionSeq = genRef.GenContentsDict[ChrNo].ToString((int) Start, (int) (End - Start)).ToCharArray();
        if (SNVs != null)
        {
            foreach (var snv in SNVs)
            {
                long loc = snv.Pos - Start;
                regionSeq[(int)loc] = snv.Alt.ToString()[0];
            }
        }
        regionSeq = Forward ? regionSeq : regionSeq.Reverse().ToArray();
        return new string(regionSeq);
    }
}