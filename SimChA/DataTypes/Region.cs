// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

// A region is zero indexed start-inclusive, end-exclusive, e.g. [0, 1) is a region of length 1 containing the first base.
public record Region(long Start, long End, ChrID ChrID, 
    Dictionary<long, SNV>? SNVDict = null, bool Forward = true) : GenRange(Start, End, ChrID.ChrNo)
{
    public static string DirToStr(bool dir) => dir ? ">" : "<";

    private string DirString => DirToStr(Forward);
    
    public override string ToString() => $"{ChrID}{DirString}[{Start}:{End})";
    public Dictionary<long, SNV> GetSNVs() => SNVDict;
}