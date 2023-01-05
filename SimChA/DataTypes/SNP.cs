namespace SimChA.DataTypes;

public struct SNP
{
    public ChrNo Chrom;
    public int Pos;
    public long AbsPos;
    public Nucleotides Ref;
    public Nucleotides Alt;
    public bool Heterozygous;
    public int Id;

    public SNP(ChrNo chrom, int pos, bool heterozygous, int id, Nucleotides re, Nucleotides alt)
    {
        Chrom = chrom;
        Pos = pos;
        AbsPos = ReferenceGenome.ChromosomeStartMap[chrom] + pos;
        Ref = re;
        Alt = alt;
        Heterozygous = heterozygous;
        Id = id;
    }

    public SNP(ChrNo chrom, int pos, bool heterozygous, int id)
    {
        Chrom = chrom;
        Pos = pos;
        AbsPos = ReferenceGenome.ChromosomeStartMap[chrom] + pos;
        Ref = Nucleotides.A;
        Alt = Nucleotides.G;
        Heterozygous = heterozygous;
        Id = id;
    }

    public override string ToString()
        => $"{Chrom} {Pos}: {Ref} / {Alt} ({(Heterozygous ? "het" : "hom")})";
}