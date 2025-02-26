using SimChA.Computation;

namespace SimChA.Data;

public class Contig
{
    private List<Region> _regions;

    public Contig()
    {
        _regions = new List<Region>();
    }

    public Contig(IEnumerable<Region> regions)
    {
        _regions = regions.Where(r => r.Length > 0).ToList();
    }

    public Contig(Contig other)
    {
        _regions = RegionOps.Copy(other._regions);
    }

    public static Contig Concat(IEnumerable<Contig> contigs)
        => new(contigs.SelectMany(c => c._regions));

    public long Length()
        => Length(_regions);
    
    public int CountRegions()
        => _regions.Count;

    public (int CNA, int CNB, int SNV) GetCNs(GenRange segRegion)
    {
        int cna = 0;
        int cnb = 0;
        int snvs = 0;
        foreach (var reg in _regions.Where(segRegion.IsInsideOf))
        {
            if (reg.Hap1)
            {
                cna += 1;
            }
            else
            {
                cnb += 1;
            }
            snvs += reg.CountSNVs(reg.Start, reg.End);
        }
        return (cna, cnb, snvs);
    }

    public Dictionary<string, List<int>> CalcBreaks()
    {
        Dictionary<string, List<int>> breaks = new();
        foreach (var region in _regions)
        {
            if (!breaks.ContainsKey(region.Chrom))
            {
                breaks[region.Chrom] = new List<int>();
            }
            breaks[region.Chrom].Add((int)region.AbsStart);
            breaks[region.Chrom].Add((int)region.AbsEnd);
        }
        return breaks;
    }

    public bool Any()
        => Length() > 0;
    
    public IEnumerable<string> GetSeq(GenRef genRef)
        => _regions.Select(r => r.GetSeq(genRef));
    
    public List<SNV> GetSNVs()
        => _regions.SelectMany(r => r.SNVs).ToList();

    public static long Length(IEnumerable<Region> regions)
        => regions.Sum(r => r.Length);

    public static string ToString(IEnumerable<Region> regions)
        => "[" + string.Join("~", regions.Select(r => r.ToString())) + "]";

    public override string ToString()
        => ToString(_regions);

    public IEnumerable<Region> FindChrRegions(string chrNo)
        => _regions.Where(r => r.Chrom == chrNo);

    public void Clear()
    {
        _regions.Clear();
    }

    public void DeleteRange(long start, long end)
    { 
        _regions = RegionOps.DeleteRange(_regions, start, end);
    }

    public Contig Split(long pos, bool keepFirst)
    {
        var (first, second) = RegionOps.SplitRegions(_regions, pos);
        _regions = keepFirst ? first : second;
        return new Contig(keepFirst ? second : first);
    }

    public void Join(Contig other)
    {
        _regions = RegionOps.ConcatRegions(_regions, other._regions);
    }

    public void InvertRange(long invStart, long invEnd)
    {
        var inverse = RegionOps.CopyRange(_regions, invStart, invEnd);
        RegionOps.Revert(inverse);
        var deleted = RegionOps.DeleteRange(_regions, invStart, invEnd);
        var (first, second) = RegionOps.SplitRegions(deleted, invStart);
        _regions = RegionOps.ConcatRegions(new[] { first, inverse, second });
    }

    public void Revert()
    {
        RegionOps.Revert(_regions);
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
            var inverse = RegionOps.Copy(second);
            RegionOps.Revert(inverse);
            _regions = RegionOps.ConcatRegions(inverse, second);
        }
        else
        {
            var inverse = RegionOps.Copy(first);
            RegionOps.Revert(inverse);
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
            RegionOps.Revert(ranges);
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
        RegionOps.PointMutateRegion(_regions, location, newNucleotide);
    }
    
    public List<string> GetPresentGenes(string chrom, List<Gene> geneList)
    {
        List<string> presentGenes = new();
        foreach (var reg in _regions.Where(r => r.Chrom == chrom && r.Forward))
        {
            int geneIndex = 0;
            while (geneIndex < geneList.Count && reg.Start > geneList[geneIndex].Start)
            {
                geneIndex++;
            }
            while (geneIndex < geneList.Count && geneList[geneIndex].End <= reg.End)
            {
                presentGenes.Add(geneList[geneIndex].Name);
                geneIndex++;
            }
        }
        return presentGenes;
    }

    public List<(long start, long end)> GetCentromeres(Dictionary<string, GenRange> centMap)
    {
        List<(long start, long end)> centromereList = new();
        long currentPos = 0;
        foreach (var reg in _regions)
        {
            var cent = centMap[reg.Chrom];
            if (cent.IsInsideOf(reg))
            {
                centromereList.Add((currentPos + cent.Start, currentPos + cent.End));
            }
            currentPos += reg.Length;
        }
        return centromereList;
    }

    public void MergeRegions()
    {
        _regions = RegionOps.MergeRegions(_regions);
    }
}
