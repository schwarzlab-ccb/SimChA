// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using System.Collections.Immutable;
using SimChA.Computation;

namespace SimChA.Data;

public class Contig
{
    private List<Region> _regions;

    public Contig()
        => _regions = new List<Region>();

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

    public List<Region> GetRegions()
        => _regions;

    public static long Length(IEnumerable<Region> regions)
        => regions.Sum(r => r.Length);

    public static string ToString(IEnumerable<Region> regions)
        => "[" + string.Join("~", regions.Select(r => r.ToString())) + "]";

    public override string ToString()
        => ToString(_regions);

    public IEnumerable<Region> FindRegionsOfChr(string chrNo)
        => _regions.Where(r => r.ChrNo == chrNo);

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

    public Contig GetSubContig(long start, long end, bool inverse = false)
    {
        var ranges = RegionOps.CopyRange(_regions, start, end);
        if (inverse)
        {
            ranges = RegionOps.InvertRegions(ranges);
        }
        return new Contig(ranges);
    }

    public void InsertContig(Contig other, long location)
    {
        if (location >= other.Length() - 1)
        {
            AppendContig(other);
        }
        else
        {
            var (first, second) = RegionOps.SplitRegions(_regions, location);
            _regions = RegionOps.ConcatRegions(new[] { first, other._regions, second });
        }
    }

    public void AppendContig(Contig other)
    {
        _regions = RegionOps.ConcatRegions(_regions, other._regions);
    }
    
    public void PointMutate(long location, Nucleotide newNucleotide)
    {
        _regions = RegionOps.PointMutateRegion(_regions, location, newNucleotide);
    }
    
    public List<Gene> GetPresentGenes(Dictionary<string, List<Gene>> geneLists)
    {
        List<Gene> presentGenes = new();
        foreach (var reg in _regions)
        {
            var geneList = geneLists[reg.ChrNo];
            if (reg.Forward && geneList.Count > 0)
            {
                int geneIndex = 0;
                while (geneIndex < geneList.Count && reg.Start > geneList[geneIndex].Range.Start)
                {
                    geneIndex++;
                }
                while (geneIndex < geneList.Count && geneList[geneIndex].Range.End <= reg.End)
                {
                    presentGenes.Add(geneList[geneIndex]);
                    geneIndex++;
                }
            }
        }
        return presentGenes;
    }

    public List<(long start, long end)> GetCentromeres(IImmutableDictionary<string, (long start, long end)> centromeres)
    {
        List<(long start, long end)> centromereList = new();
        long currentPos = 0;
        foreach (var reg in _regions)
        {
            (long centStart, long centEnd) = centromeres[reg.ChrNo];
            if (reg.Start <= centStart && reg.End >= centEnd)
            {
                centromereList.Add((currentPos + centStart, currentPos + centEnd));
            }
            currentPos += reg.Length;
        }
        return centromereList;
    }

    public void MergeRegions()
        => _regions = RegionOps.MergeRegions(_regions);
}
