namespace SimChA.DataTypes;

public record SNVRegion(long Start, long End, ChrNo ChrNo, bool Hap1, Dictionary<long, SNV> SNVDict, bool Forward = true) 
    : Region(Start, End, ChrNo, Hap1, Forward)
{
    // TODO: @cody should add SNV output
    public override string ToString() => base.ToString();
}