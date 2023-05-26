namespace SimChA.DataTypes;

public record GeneRange(long Start, long End, ChrNo ChrNo)
{
    public long Length => End - Start;
    
    public bool IsInside(GeneRange other) 
        => Start >= other.Start && End <= other.End && ChrNo == other.ChrNo;

    public bool Overlaps(GeneRange other)
        => Start < other.End && End > other.Start && ChrNo == other.ChrNo;
    
    public override string ToString() 
        => $"{ChrNo}[{Start}:{End})";
        
    public bool IsInside(long start, long end, ChrNo chrNo)
        => Start >= start && End <= end && ChrNo == chrNo;
}