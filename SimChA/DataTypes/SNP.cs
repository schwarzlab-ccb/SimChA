namespace SimChA.DataTypes;

public record SNP
{
    private readonly Nucleotides _ref;
    private readonly Nucleotides _alt;
    public readonly ChrNo ChrNo;
    public readonly int Pos;
    public readonly long AbsPos;
    public readonly bool Heterozygous;
    public readonly int Id;

    public SNP(ChrNo chrNo, int pos, bool heterozygous, int id)
    {
        ChrNo = chrNo;
        Pos = pos;
        AbsPos = ReferenceGenome.ChromosomeStartMap[chrNo] + pos;
        _ref = Nucleotides.A;
        _alt = Nucleotides.G;
        Heterozygous = heterozygous;
        Id = id;
    }

    public override string ToString()
        => $"{ChrNo} {Pos}: {_ref} / {_alt} ({(Heterozygous ? "het" : "hom")})";
}