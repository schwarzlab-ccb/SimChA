namespace SimChA.Data;

public class CopyNumber : GenRange
{
    public int CNH1 { get; }
    public int CNH2 { get; }
    public int NSNVs { get; }
    
    public CopyNumber(long start, long end, string chrom, int cnh1, int cnh2, int nsnvs)
        : base(start, end, chrom)
    {
        CNH1 = cnh1;
        CNH2 = cnh2;
        NSNVs = nsnvs;
    }
    
    private static string NaIfNegStr(int val)
        => val < 0 ? "NA" : $"{val}";

    public static string Header()
        => "chrom\tstart\tend\tcn_a\tcn_b\tn_snvs";
    
    public string ToTSV()
        => string.Join('\t', Chrom, Start + 1, End, NaIfNegStr(CNH1), NaIfNegStr(CNH2), NaIfNegStr(NSNVs));
}
