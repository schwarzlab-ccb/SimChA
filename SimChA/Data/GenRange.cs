namespace SimChA.Data;

public record GenRange(long Start, long End, string Chrom)
{
    public long Length => End - Start;
    
    public bool Forward => Start >= 0;

    public long AbsStart => Math.Max(Start, -End);
    
    public long AbsEnd => Math.Max(End, -Start);

    public static string DirToStr(bool dir) => dir ? ">" : "<";

    private string DirString => DirToStr(Forward);
    
    // True if this range is inside the other range
    public bool IsInsideOf(GenRange other) 
        => AbsStart >= other.AbsStart && AbsEnd <= other.AbsEnd && Chrom == other.Chrom;

    // True if this range shares at least one position with the other range
    public bool Overlaps(GenRange other)
        => Start < other.End && End > other.Start && Chrom == other.Chrom;
    
    public override string ToString() 
        => $"{DirString}{Chrom}[{Start}:{End})";
}
