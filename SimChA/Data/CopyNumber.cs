namespace SimChA.Data;

public record CopyNumber(long Start, long End, string Chrom, int CNH1, int CNH2, int NSNVs) 
    : GenRange(Start, End, Chrom)
{
    public CopyNumber(GenRange range, int cnH1, int cnH2, int nSNVs) 
        : this(range.Start, range.End, range.Chrom, cnH1, cnH2, nSNVs) { }
    
    private static string NAIfNegStr(int val)
        => val < 0 ? "NA" : $"{val}";
    
    public string ToTSV()
        => string.Join('\t', Chrom, Start + 1, End, NAIfNegStr(CNH1), NAIfNegStr(CNH2), NAIfNegStr(NSNVs));
}
