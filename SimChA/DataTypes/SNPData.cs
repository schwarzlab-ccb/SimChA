namespace SimChA.DataTypes;

public struct SNPData
{
    public SNP Snp;
    public float Logr;
    public float Baf;

    public SNPData(SNP snp, float logr, float baf)
    {
        Snp = snp;
        Logr = logr;
        Baf = baf;
    }

    public string PrintBAF()
        => $"{Snp.Id}\t{Snp.Chrom}\t{Snp.Pos}\t{Baf}";

    public string PrintLogR()
        => $"{Snp.Id}\t{Snp.Chrom}\t{Snp.Pos}\t{Logr}";
}