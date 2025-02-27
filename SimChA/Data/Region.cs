namespace SimChA.Data;

// A region is zero indexed start-inclusive, end-exclusive, e.g. [0, 1) is a region of length 1 containing the first base.
public class Region : GenRange
{
    public bool Hap1 { get; }
    public List<SNV> SNVs { get; }
    
    public Region(long start, long end, string chrom, bool hap1, List<SNV> snvs) : base(start, end, chrom)
    {
        Hap1 = hap1;
        SNVs = snvs;
    }

    public Region(Region other) : base(other)
    {
        Hap1 = other.Hap1;
        SNVs = new List<SNV>(other.SNVs);
    }

    public override bool Equals(object? obj) 
        => obj is Region other && base.Equals(other) && Hap1 == other.Hap1 && SNVs.SequenceEqual(other.SNVs);
    
    public void ResizeFront(long howMuch)
    {
        Start += howMuch;
        UpdateSNVs();
    }
    
    public void ResizeBack(long howMuch)
    {
        End = Start + howMuch;
        UpdateSNVs();
    }    
    
    public void ResizeBoth(long start, long end)
    {
        End = Start + end;
        Start += start;
        UpdateSNVs();
    }
    
    private void UpdateSNVs()
    {
        foreach (var snv in SNVs.Where(snv => snv.Pos <= AbsStart || AbsEnd <= snv.Pos).ToList())
        {
            SNVs.Remove(snv);
        }
    }

    private void AddSNVs(List<SNV> snvs)
    {
        SNVs.AddRange(snvs);
    }

    public void AddSNV(long offset, Nucleotide newNucleotide)
    {
        int index = SNVs.FindIndex(s => s.Pos == AbsStart + offset);
        if (index >= 0)
        {
            // Update the existing SNV
            SNVs[index] = SNVs[index] with { Alt = newNucleotide };
        }
        else
        {
            SNVs.Add(new SNV(AbsStart + offset, Chrom, newNucleotide));
        }
    }
    
    private static string HapToString(bool parent) 
        => parent ? "H1" : "H2";

    private string HapString 
        => HapToString(Hap1);
    
    public override string ToString() 
        => $"{HapString}{base.ToString()}";

    public int CountSNVs(long start, long end)
        => SNVs?.Count(s => s.Pos >= start && s.Pos < end) ?? 0;

    public string GetSeq(GenRef genRef)
    {
        char[] regionSeq = genRef.GenContentsDict[Chrom].ToString((int) Start, (int) (End - Start)).ToCharArray();
        if (SNVs != null)
        {
            foreach (var snv in SNVs)
            {
                long loc = snv.Pos - Start;
                regionSeq[(int)loc] = snv.Alt.ToString()[0];
            }
        }
        regionSeq = Forward ? regionSeq : regionSeq.Reverse().ToArray();
        return new string(regionSeq);
    }

    public void MergeWith(Region cur)
    {
        End = cur.End;
        AddSNVs(cur.SNVs);
    }
}