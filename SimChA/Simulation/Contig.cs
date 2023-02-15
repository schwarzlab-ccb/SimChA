// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.Misc;

namespace SimChA.Simulation;

public class Contig
{
    private List<Region> _regions;

    public Contig(Region initialRegion) 
        => _regions = new List<Region> { initialRegion };

    public Contig(IEnumerable<Region> regions)
        => _regions = regions.Where(r => r.Length > 0).ToList();

    public Contig(Contig other) 
        => _regions = new List<Region>(other._regions);

    public long Length() 
        => Length(_regions);

    public bool Any() 
        => Length() > 0;

    public static long Length(IEnumerable<Region> regions)
        => regions.Sum(r => r.Length);

    public static string ToString(IEnumerable<Region> regions)
        => "[" + string.Join(",", regions.Select(r => r.ToString())) + "]";

    public override string ToString()
        => ToString(_regions);
    
    public IEnumerable<Region> FindRegionsOfChr(ChrNo chrNo)
        => _regions.Where(r => r.ChrID.ChrNo == chrNo);

    public void Clear()
        => _regions.Clear();

    public void DeleteRange(long start, long end)
        => _regions = RegionOps.DeleteRange(_regions, start, end);

    public Contig Split(long pos, bool keepFirst)
    {
        var (first, second) = RegionOps.SplitRegions(_regions, pos);
        _regions = keepFirst ? first : second;
        return new Contig(keepFirst ? second : first);
    }

    public void Join(Contig other)
        => _regions = RegionOps.ConcatRegions(_regions, other._regions);

    public void InvertRange(long invStart, long invEnd)
    {
        var copy = RegionOps.CopyRange(_regions, invStart, invEnd);
        var inverse = RegionOps.InvertRegions(copy);
        var deleted = RegionOps.DeleteRange(_regions, invStart, invEnd);
        var (first, second) = RegionOps.SplitRegions(deleted, invStart);
        _regions = RegionOps.ConcatRegions(new[] { first, inverse, second });
    }

    public void DuplicateRange(long start, long end)
    {
        var copy = RegionOps.CopyRange(_regions, start, end);
        var (first, second) = RegionOps.SplitRegions(_regions, start);
        _regions = RegionOps.ConcatRegions(new[] { first, copy, second });
    }

    public void Bridge(long pos, bool cutFront)
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

    public void ScatterAndGather(List<long> locs, IEnumerable<int> indices)
    {
        // First region
        var newRegions = new List<List<Region>> { RegionOps.CopyRange(_regions, 0, locs[0]) };
        // Internal regions
        for (int i = 0; i < locs.Count - 1; i++)
        {
            long start = locs[i];
            long end = locs[i + 1];
            var copy = RegionOps.CopyRange(_regions, start, end);
            newRegions.Add(copy);
        }
        // Last region
        newRegions.Add(RegionOps.CopyRange(_regions, locs.Last(), Length(_regions)));
        
        var selectedRegions = indices.Select(i => newRegions[i]);
        _regions = RegionOps.ConcatRegions(selectedRegions);
    }
    
    public IEnumerable<Gene> GetPresentGenes(Dictionary<ChrNo, List<Gene>> geneLists)
        => _regions.SelectMany(r => geneLists[r.ChrID.ChrNo].FindAll(g => g.Region.IsInside(r)));
}