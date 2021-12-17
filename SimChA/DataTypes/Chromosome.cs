// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.Simulation;

namespace SimChA.DataTypes;

public class Chromosome
{
    private List<Region> _regions;

    public Chromosome(Region initialRegion)
    {
        _regions = new List<Region> {initialRegion};
    }

    public Chromosome(IEnumerable<Region> regions)
    {
        _regions = regions.Where(r => r.Length > 0).ToList();
    }

    public Chromosome(Chromosome other)
    {
        _regions = new List<Region>(other._regions);
    }

    public int Length() => Length(_regions);

    public static int Length(List<Region> regions) => regions.Sum(r => r.Length);


    public static string ToString(List<Region> regions)
        => "[" + string.Join(",", regions.Select(r => r.ToString())) + "]";

    public override string ToString() => ToString(_regions);

    public void DeleteRange(int start, int end)
    {
        _regions = ChrMutations.DeleteRange(_regions, start, end);
    }

    public void DuplicateRange(int start, int end)
    {
        var newRegions = ChrMutations.CopyRange(_regions, start, end);
    }
}