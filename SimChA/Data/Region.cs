namespace SimChA.Data;

// A region is zero indexed start-inclusive, end-exclusive, e.g. [0, 1) is a region of length 1 containing the first base.
public class Region : GenRange
{
    public bool Hap1 { get; }
    public List<SNV> SNVs { get; }
    public List<Gene> Genes { get; }
    
    public Region(long start, long end, string chrom, bool hap1, 
        List<SNV> snvs, List<Gene> presentGenes) : base(start, end, chrom)
    {
        Hap1 = hap1;
        SNVs = snvs;
        Genes = presentGenes;
    }

    public Region(Region other) : base(other)
    {
        Hap1 = other.Hap1;
        SNVs = [..other.SNVs];
        Genes = [..other.Genes];
    }

    public override bool Equals(object? obj) 
        => obj is Region other && base.Equals(other) && Hap1 == other.Hap1 && SNVs.SequenceEqual(other.SNVs);
    
    public void ResizeFront(long howMuch)
    {
        Start += howMuch;
        UpdateRegion();
    }
    
    public void ResizeBack(long howMuch)
    {
        End = Start + howMuch;
        UpdateRegion();
    }    
    
    public void ResizeBoth(long start, long end)
    {
        End = Start + end;
        Start += start;
        UpdateRegion();
    }

    private void UpdateRegion()
    {
        SNVs.RemoveAll(snv => snv.Pos <= AbsStart || AbsEnd <= snv.Pos);
        Genes.RemoveAll(g => !g.IsInsideOf(this));
    }
    
    public void AddSNV(long offset, Nucleotide oldNucleotide, Nucleotide newNucleotide)
    {
        int index = SNVs.FindIndex(s => s.Pos == AbsStart + offset);
        if (index >= 0)
        {
            // Update the existing SNV if newNucleotide is different from ref
            if (newNucleotide != oldNucleotide)
            {
                SNVs[index] = SNVs[index] with { Alt = newNucleotide };
            }
            else
            {
                SNVs.RemoveAt(index);
            }
            
        }
        else
        {
            SNVs.Add(new SNV(AbsStart + offset, Chrom, oldNucleotide, newNucleotide));
        }
    }
    
    private static string HapToString(bool parent) 
        => parent ? "H1" : "H2";

    private string HapString 
        => HapToString(Hap1);
    
    public override string ToString() 
        => $"{HapString}{base.ToString()}";

    public int CountSNVs(long start, long end)
        => SNVs.Count(s => s.Pos >= start && s.Pos < end);

    public string GetSeq(GenRef genRef)
    {
        char[] regionSeq = genRef.GetGenContents(Chrom, (int) Start, (int) (End-Start));
        // If the GenContents haven't been set, we return the null case
        if (regionSeq is ['N']) 
        {
            return new string(regionSeq);
        }
        foreach (var snv in SNVs)
        {
            long loc = snv.Pos - Start;
            regionSeq[(int)loc] = snv.Alt.ToString()[0];
        }
        regionSeq = Forward ? regionSeq : regionSeq.Reverse().ToArray();
        return new string(regionSeq);
    }

    public void MergeWithNext(Region next)
    {
        End = next.End;
        SNVs.AddRange(next.SNVs);
        Genes.AddRange(next.Genes);
    }
}