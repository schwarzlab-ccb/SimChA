// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.DataTypes;
namespace SimChA.Simulation;

public class Chromosome
{
    private List<Region> _regions;

    public Chromosome(Region initialRegion)
    {
        _regions = new List<Region> { initialRegion };
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

    public static int Length(List<Region> regions) 
        => regions.Sum(r => r.Length);

    public static string ToString(List<Region> regions)
        => "[" + string.Join(",", regions.Select(r => r.ToString())) + "]";

    public override string ToString() 
        => ToString(_regions);

    public List<Region> GetAllRegions() 
        => _regions.ToList();

    public List<Region> GetRegions(ChromID chromID) 
        => _regions.Where(r => r.ChromId.Equals(chromID)).ToList();

    public void DeleteRange(int start, int end)
    {
        _regions = RegionOps.DeleteRange(_regions, start, end);
    }

    public Chromosome Split(int pos, bool keepFirst)
    {
        var (first, second) = RegionOps.SplitRegions(_regions, pos);
        _regions = keepFirst ? first : second;
        return new Chromosome(keepFirst ? second : first);
    }

    public void Join(Chromosome other, bool prepend)
    {
        if (prepend)
        {
            _regions = RegionOps.ConcatRegions(other._regions, _regions);
        }
    }

    public void InvertRange(int invStart, int invEnd)
    {
        var copy = RegionOps.CopyRange(_regions, invStart, invEnd);
        var inverse = RegionOps.InvertRegions(copy);
        var deleted = RegionOps.DeleteRange(_regions, invStart, invEnd);
        var (first, second) = RegionOps.SplitRegions(deleted, invStart);
        _regions = RegionOps.ConcatRegions(new[] { first, inverse, second });
    }
    
    public void DuplicateRange(int start, int end)
    {
        var copy = RegionOps.CopyRange(_regions, start, end);
        var (first, second) = RegionOps.SplitRegions(_regions, start);
        _regions = RegionOps.ConcatRegions(new[] { first, copy, second });
    }

    public void Bridge(int pos, bool cutFront)
    {
        var (first, second) = RegionOps.SplitRegions(_regions, pos);
        if (cutFront)
        {
            var inverse = RegionOps.InvertRegions(second);
            _regions = RegionOps.ConcatRegions(inverse, second);
        }
        else
        {
            var inverse = RegionOps.InvertRegions(first);
            _regions = RegionOps.ConcatRegions(first, inverse);
        }
    }

    public void ScatterAndGather(List<int> locs, int count)
    {
        var newRegions = new List<List<Region>> { RegionOps.CopyRange(_regions, 0, locs[0]) };
        for (int i = 0; i < locs.Count - 1; i++)
        {
            int start = locs[i];
            int end = locs[i + 1];
            var copy = RegionOps.CopyRange(_regions, start, end);
            newRegions.Add(copy);
        }
        newRegions.Add(RegionOps.CopyRange(_regions, locs.Last(), Length(_regions)));
        
        var keptRegions = newRegions.Shuffle().Take(count);
        _regions = RegionOps.ConcatRegions(keptRegions);
    }
}