namespace SimChA.Data;

public record GenRange(long Start, long End, string ChrNo)
{
    public long Length => End - Start;
    
    public bool Forward => Start >= 0;
    
    // True if this range is inside the other range
    public bool IsInside(GenRange other) 
        => Math.Max(Start, -End) >= other.Start && Math.Max(End, -Start) <= other.End && ChrNo == other.ChrNo;

    // True if this range shares at least one position with the other range
    public bool Overlaps(GenRange other)
        => Start < other.End && End > other.Start && ChrNo == other.ChrNo;
    
    public override string ToString() 
        => $"{ChrNo}[{Start}:{End})";
}
