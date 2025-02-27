using SimChA.Computation;

namespace SimChA.Data;

public class Contig
{
    private long _length;
    private List<Region> _regions;
    private List<Region> Regions
    {
        get => _regions;
        set
        {
            _regions = value;
            _length = _regions.Count > 0 ? RegionOps.CountLength(_regions) : 0;
        }
    }

    public Contig()
    {
        _regions = new List<Region>();
        _length = 0;
    }

    public Contig(IEnumerable<Region> regions)
    {
        _regions = regions.Where(r => r.Length > 0).ToList();
        _length = RegionOps.CountLength(_regions);
    }

    public Contig(Contig other)
    {
        _regions = RegionOps.Copy(other.Regions);
        _length = other.Length;
    }

    public static Contig Concat(IEnumerable<Contig> contigs)
        => new(contigs.SelectMany(c => c.Regions));
    
    public long Length
        => _length;
    
    public int CountRegions()
        => Regions.Count;

    public (int CNA, int CNB, int SNV) GetCNs(GenRange segRegion)
    {
        int cna = 0;
        int cnb = 0;
        int snvs = 0;
        foreach (var reg in Regions.Where(segRegion.IsInsideOf))
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
        foreach (var region in Regions)
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
        => _length > 0;
    
    public IEnumerable<string> GetSeq(GenRef genRef)
        => Regions.Select(r => r.GetSeq(genRef));
    
    public List<SNV> GetSNVs()
        => Regions.SelectMany(r => r.SNVs).ToList();

    public static string ToString(IEnumerable<Region> regions)
        => "[" + string.Join("~", regions.Select(r => r.ToString())) + "]";

    public override string ToString()
        => ToString(Regions);

    public IEnumerable<Region> FindChrRegions(string chrNo)
        => Regions.Where(r => r.Chrom == chrNo);

    public void Clear()
    {
        Regions.Clear();
        _length = 0;
    }

    public void DeleteRange(long start, long end)
    { 
        Regions = RegionOps.DeleteRange(Regions, start, end);
    }

    public Contig Split(long pos, bool keepFirst)
    {
        var (first, second) = RegionOps.SplitRegions(Regions, pos);
        Regions = keepFirst ? first : second;
        return new Contig(keepFirst ? second : first);
    }

    public void Join(Contig other)
    {
        Regions = RegionOps.ConcatRegions(Regions, other.Regions);
    }

    public void InvertRange(long invStart, long invEnd)
    {
        var inverse = RegionOps.CopyRange(Regions, invStart, invEnd);
        RegionOps.Revert(inverse);
        var deleted = RegionOps.DeleteRange(Regions, invStart, invEnd);
        var (first, second) = RegionOps.SplitRegions(deleted, invStart);
        Regions = RegionOps.ConcatRegions(new[] { first, inverse, second });
    }

    public void Revert()
    {
        RegionOps.Revert(Regions);
    }

    public void DuplicateRange(long start, long end)
    {
        var copy = RegionOps.CopyRange(Regions, start, end);
        var (first, second) = RegionOps.SplitRegions(Regions, start);
        Regions = RegionOps.ConcatRegions(new[] { first, copy, second });
    }

    public void Bridge(long pos, bool cutFront)
    {
        var (first, second) = RegionOps.SplitRegions(Regions, pos);
        if (cutFront)
        {
            var inverse = RegionOps.Copy(second);
            RegionOps.Revert(inverse);
            Regions = RegionOps.ConcatRegions(inverse, second);
        }
        else
        {
            var inverse = RegionOps.Copy(first);
            RegionOps.Revert(inverse);
            Regions = RegionOps.ConcatRegions(first, inverse);
        }
    }

    public IEnumerable<Contig> Scatter(List<long> locs)
    {
        var regions = RegionOps.Scatter(locs, Regions);
        return regions.Select(r => new Contig(r));
    }

    public void ScatterAndGather(List<long> locs, IEnumerable<int> indices)
    {
        var newRegions = RegionOps.Scatter(locs, Regions);
        Regions = RegionOps.Gather(newRegions, indices);
    }

    public Contig GetSubContig(long start, long end, bool inverse = false)
    {
        var ranges = RegionOps.CopyRange(Regions, start, end);
        if (inverse)
        {
            RegionOps.Revert(ranges);
        }
        return new Contig(ranges);
    }

    public void InsertContig(Contig other, long location)
    {
        if (location >= other.Length - 1)
        {
            AppendContig(other);
        }
        else
        {
            var (first, second) = RegionOps.SplitRegions(Regions, location);
            Regions = RegionOps.ConcatRegions(new[] { first, other.Regions, second });
        }
    }

    public void AppendContig(Contig other)
    {
        Regions = RegionOps.ConcatRegions(Regions, other.Regions);
    }
    
    public void PointMutate(long location, Nucleotide newNucleotide)
    {
        RegionOps.PointMutateRegion(Regions, location, newNucleotide);
    }
    
    public List<string> GetPresentGenes(string chrom, List<Gene> geneList)
    {
        List<string> presentGenes = new();
        foreach (var reg in Regions.Where(r => r.Chrom == chrom && r.Forward))
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
        foreach (var reg in Regions)
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
        Regions = RegionOps.MergeRegions(Regions);
    }
}
