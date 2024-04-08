// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

// q-arm of a chromosome is just a dedicated region
public record QArm(long Start, long End, string ChrNo, bool Hap1, bool Forward = true, Dictionary<long, Nucleotide>? SNVDict = null) : Region(Start, End, ChrNo, Hap1, Forward, SNVDict)
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
