namespace SimChA.DataTypes;

public struct SNPData
{
    public SNP Snp;
    public float LogR;
    public float Baf;

    public SNPData(SNP snp, float logR, float baf)
    {
        Snp = snp;
        LogR = logR;
        Baf = baf;
    }

    public string PrintBAF()
        => $"{Snp.Id}\t{Snp.ChrNo}\t{Snp.Pos}\t{Baf}";

    public string PrintLogR()
        => $"{Snp.Id}\t{Snp.ChrNo}\t{Snp.Pos}\t{LogR}";
}