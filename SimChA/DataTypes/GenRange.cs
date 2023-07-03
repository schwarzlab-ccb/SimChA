namespace SimChA.DataTypes;

public record GenRange(long Start, long End, string ChrNo)
{
    public long Length => End - Start;
    
    // True if this range is inside the other range
    public bool IsInside(GenRange other) 
        => Start >= other.Start && End <= other.End && ChrNo == other.ChrNo;

    // True if this range shares at least one position with the other range
    public bool Overlaps(GenRange other)
        => Start < other.End && End > other.Start && ChrNo == other.ChrNo;
    
    public override string ToString() 
        => $"{ChrNo}[{Start}:{End})";
}