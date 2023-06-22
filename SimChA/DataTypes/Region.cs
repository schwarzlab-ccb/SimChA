// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

// A region is zero indexed start-inclusive, end-exclusive, e.g. [0, 1) is a region of length 1 containing the first base.
// TODO: The direction should not be part of the region, should be part of the chromosome
public record Region(long Start, long End, ChrNo ChrNo, bool Hap1, bool Forward = true, Dictionary<long, Nucleotide>? SNVDict = null) : GenRange(Start, End, ChrNo)
{ 
    private static string HapToString(bool parent) 
        => parent ? "H1" : "H2";
    
    private string HapString => HapToString(Hap1);

    public static string DirToStr(bool dir) => dir ? ">" : "<";

    private string DirString => DirToStr(Forward);
    
    public override string ToString() => $"{HapString}{DirString}{ChrNo}[{Start}:{End})";

    public int NumSNVsBetween(long start, long end)
        => (SNVDict!=null) ? SNVDict.Keys.Count(loc => loc>= start && loc <= end) : 0 ;
}
