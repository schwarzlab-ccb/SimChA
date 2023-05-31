namespace SimChA.DataTypes;

using SimChA.Computation;
public record GeneRegionList : GeneRange
{
    public Dictionary<GeneListType, List<Gene>> GeneLists { get; init; }

    public GeneRegionList(long start, long end, ChrNo chrNo) : base(start, end, chrNo)
    {
        GeneLists = Fitness.GetGeneList(start, end, chrNo);
    }
}