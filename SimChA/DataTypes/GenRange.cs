namespace SimChA.DataTypes;

public record GenRange(long Start, long End, ChrNo ChrNo)
{
    public long Length => End - Start;
    
    public bool IsInside(GenRange other) 
        => Start >= other.Start && End <= other.End && ChrNo == other.ChrNo;

    public bool Overlaps(GenRange other)
        => Start < other.End && End > other.Start && ChrNo == other.ChrNo;
    
    public override string ToString() 
        => $"{ChrNo}[{Start}:{End})";
}