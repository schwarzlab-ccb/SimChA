// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

// A region is zero indexed start-inclusive, end-exclusive, e.g. [0, 1) is a region of length 1 containing the first base.
// TODO: The direction should not be part of the region, should be part of the chromosome
public record Region(long Start, long End, ChrID ChrID, bool Forward = true) : GeneRegionList(Start, End, ChrID.ChrNo)
{
    private string DirString => Forward ? ">" : "<";
    
    public override string ToString() => $"{ChrID}{DirString}[{Start}:{End})";
}