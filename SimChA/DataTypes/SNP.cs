namespace SimChA.DataTypes;

public record SNP(ChrNo ChrNo, int Pos, bool Heterozygous, int Id)
{
    public readonly ChrNo ChrNo = ChrNo;
    public readonly int Pos = Pos;
    public long AbsPos => ReferenceGenome.ChromosomeStartMap[ChrNo] + Pos;
    public readonly bool Heterozygous = Heterozygous;
    public readonly int Id = Id; // Id is the index in the list of SNPs

    public override string ToString()
        => $"{ChrNo} {Pos}: ({(Heterozygous ? "het" : "hom")})";
}