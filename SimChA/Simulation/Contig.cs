// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.Computation;
using SimChA.DataTypes;

namespace SimChA.Simulation;

public class Contig
{
    public bool Circular { get; set; }
    
    private List<Region> _regions;

    public Contig(Region initialRegion) 
        => _regions = new List<Region> { initialRegion };

    public Contig(IEnumerable<Region> regions)
        => _regions = regions.Where(r => r.Length > 0).ToList();

    public Contig(Contig other) 
        => _regions = new List<Region>(other._regions);
    
    public static Contig Concat(IEnumerable<Contig> contigs) 
        => new(contigs.SelectMany(c => c._regions));

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

    public void Invert()
    {
        _regions = RegionOps.InvertRegions(_regions);
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

    public IEnumerable<Contig> Scatter(List<long> locs)
    {
        var regions = RegionOps.Scatter(locs, _regions);
        return regions.Select(r => new Contig(r));
    }

    public void ScatterAndGather(List<long> locs, IEnumerable<int> indices)
    {
        var newRegions = RegionOps.Scatter(locs, _regions);
        _regions = RegionOps.Gather(newRegions, indices);
    }
    
    public List<Region> GetRandomRegion(long start, long end)
        => RegionOps.CopyRange(_regions, start, end);   

    public void AddRegions(List<Region> regions, long insertion)
    {
        var (first, second) = RegionOps.SplitRegions(_regions, insertion);
        _regions = RegionOps.ConcatRegions(new[] {first, regions, second});
    }

    public void GlueNeighbours()
        => _regions = RegionOps.GlueNeighbours(_regions);

    public List<Region> GetSplitRegion(long pos, bool foreward)
    {
        var (first, second) = RegionOps.SplitRegions(_regions, pos);
        return foreward ? second : first;
    }

    public IEnumerable<Gene> GetPresentGenes(Dictionary<ChrNo, List<Gene>> geneLists)
        => _regions.SelectMany(r => geneLists[r.ChrID.ChrNo].FindAll(g => r.Forward && g.Range.IsInside(r)));
}