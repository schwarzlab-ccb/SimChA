namespace SimChA.DataTypes;

public struct SNP
{
    public ChromNum Chrom;
    public int Pos;
    public long AbsPos;
    public Nucleotides Ref;
    public Nucleotides Alt;
    public bool Heterozygous;

    public SNP(ChromNum chrom, int pos, bool heterozygous, Nucleotides re, Nucleotides alt)
    {
        Chrom = chrom;
        Pos = pos;
        AbsPos = ReferenceGenome.ChromosomeAbsoluteStart(chrom) + pos;
        Ref = re;
        Alt = alt;
        Heterozygous = heterozygous;
    }
    public SNP(ChromNum chrom, int pos, bool heterozygous)
    {
        Chrom = chrom;
        Pos = pos;
        AbsPos = ReferenceGenome.ChromosomeAbsoluteStart(chrom) + pos;
        Ref = Nucleotides.A;
        Alt = Nucleotides.G;
        Heterozygous = heterozygous;
    }

    public override string ToString()
        => $"{Chrom} {Pos}: {Ref} / {Alt} ({(Heterozygous ? "het" : "hom")})";
}