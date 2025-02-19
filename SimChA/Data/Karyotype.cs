using Extreme.Mathematics;
using SimChA.Computation;
using SimChA.IO;

namespace SimChA.Data;

public class Karyotype
{
    private GenRef GenRef { get; }
    public SexType Sex { get; }
    public double FitnessVal { get; private set; }
    // NOTE: Empty contigs are retained in the list, but not reported. This way the initial indexing is preserved.
    private readonly List<Contig> _contigs;
    
    public Karyotype(GenRef genRef, SexType sex)
    {
        GenRef = genRef;
        _contigs = genRef.GetGenotype(sex).Select(region => new Contig(region)).ToList();
        Sex = sex;
    }
    
    public Karyotype(Karyotype other)
    {
        GenRef = other.GenRef;
        _contigs = other._contigs.Select(ch => new Contig(ch)).ToList();
        Sex = other.Sex;
        FitnessVal = other.FitnessVal;
    }
    
    public Karyotype(GenRef genRef, List<Contig> contigs, SexType sex)
    {
        GenRef = genRef;
        _contigs = contigs;
        Sex = sex;
    }

    public int CountContigs() 
        => _contigs.Count(c => c.Any());
    
    public long GenomeLen()
        => _contigs.Sum(c => c.Length());
    
    public IEnumerable<SNV> GetSNVs()
        => _contigs.SelectMany(c => c.GetSNVs()).ToList();

    public Dictionary<string, List<int>> CalcBreaks()
    {
        var breakSets = GenRef.ChrIDsForSex(Sex).ToDictionary(c => c, c => new HashSet<int> {0, GenRef.ChrLengths[c] });
        foreach (var contig in _contigs)
        {
            foreach ((string chrom, var breaks) in contig.CalcBreaks())
            {
                breakSets[chrom].UnionWith(breaks);
            }
        }
        return breakSets.ToDictionary(k => k.Key, v => v.Value.OrderBy(i => i).ToList());
    }
    
    public IEnumerable<CopyNumber> CalcChrCopyNumbers(List<int> breaks, string chrom)
    {
        var result = new List<CopyNumber>();
        for (int i = 0; i < breaks.Count - 1; i++)
        {
            int start = breaks[i];
            int end = breaks[i + 1];
            var seg = new GenRange(start, end, chrom);
            var cns = _contigs.Select(c => c.GetCNs(seg));
            (int cnA, int cnB, int nSNVs) = cns.Aggregate((CNA: 0, CNB: 0, SNV: 0), (acc, vals) 
                => (acc.CNA + vals.CNA, acc.CNB + vals.CNB, acc.SNV + vals.SNV));
            var cn = new CopyNumber(seg, cnA, cnB, nSNVs);
            result.Add(cn);
        }
        return result;
    }

    public IEnumerable<int> ContigIds() 
        => _contigs.Select((c, i) => (c, i)).Where(t => t.c.Any()).Select(t => t.i);
    
    public override string ToString()
        => CountContigs() > 0 ? "[" + string.Join(";", _contigs.Where(c => c.Any())) + "]" : "[]";
    
    // @CODY this should be made private, needs to update tests
    public IEnumerable<Region> FindChrRegions(string chrNo) 
        => _contigs.SelectMany(c => c.FindChrRegions(chrNo));

    private static long GetTailSplitPos(long segLength, Contig contig, bool fiveToThree) 
        => fiveToThree ? segLength : contig.Length() - segLength;
    
    public long ContigLen(int contigId)
        => contigId < _contigs.Count ? _contigs[contigId].Length() : 0;
    
    public List<(long start, long end)> GetCentromeres(int contigId)
        => _contigs[contigId].GetCentromeres(GenRef.Centromeres);
    
    private static (long start, long end) GetIndices(Contig contig, long position, bool fiveToThree)
        => fiveToThree ? (0, position) : (position, contig.Length());
    
    // @CODY this should be made private, needs to update test
    public Contig GetContig(int contigID) 
        => _contigs[contigID];

    public IEnumerable<string> GetSeq()
        => _contigs.SelectMany(c => c.GetSeq(GenRef).Concat(new List<string> {"\n"}));
    
    public List<Gene> GetPresentGenes(Dictionary<string, List<Gene>> geneLists)
        => _contigs.SelectMany(c => c.GetPresentGenes(geneLists)).ToList();

    public double UpdateFitness(GenRef genRef, FitParams fParams)
        => FitnessVal = Fitness.Calculate(this, genRef, fParams);
    
    public void ApplyTailDeletion(int contigID, long tailLen, bool fiveToThree)
    {
        var contig = _contigs[contigID];
        long tailSplit = GetTailSplitPos(tailLen, contig, fiveToThree);
        (long tailStart, long tailEnd) = GetIndices(contig, tailSplit, fiveToThree);
        contig.DeleteRange(tailStart, tailEnd);
    }

    public void ApplyTailDuplication(int contigID, long tailLen, bool fiveToThree)
    {
        var contig = _contigs[contigID];
        long tailSplit = GetTailSplitPos(tailLen, contig, fiveToThree);
        (long tailStart, long tailEnd) = GetIndices(contig, tailSplit, fiveToThree);
        var newTail = new Contig(contig.GetSubContig(tailStart, tailEnd));
        _contigs.Add(newTail);
    }

    public void ApplyBFB(int contigID, long tailLen, bool fiveToThree)
    {
        var contig = _contigs[contigID];
        long tailSplit = GetTailSplitPos(tailLen, contig, fiveToThree);
        contig.Bridge(tailSplit, fiveToThree);
    }
    
    public void ApplyContigDeletion(int contigID)
    {
        var contig = _contigs[contigID];
        contig.Clear();
    }
    
    public void ApplyContigDuplication(int contigID)
    {
        var contig = _contigs[contigID];
        _contigs.Add(new Contig(contig));
    }

    public void ApplyInternalDuplication(int contigID, long startPos, long endPos)
    {
        var contig = _contigs[contigID];
        contig.DuplicateRange(startPos, endPos);
    }
    
    public void ApplyInvertedDuplication(int contigID, long startPos, long endPos)
    {
        var contig = _contigs[contigID];
        contig.DuplicateRange(startPos, endPos);
        contig.InvertRange(endPos, endPos + (endPos - startPos));
    }
    
    public void ApplyInternalInversion(int contigID, long startPos, long endPos)
    {
        var contig = _contigs[contigID];
        contig.InvertRange(startPos, endPos);
    }

    public void ApplyInternalDeletion(int contigID, long startPos, long endPos)
    {
        var contig = _contigs[contigID];
        contig.DeleteRange(startPos, endPos);
    }

    // Translocation might invert based on the orientation of the holiday Junction https://en.wikipedia.org/wiki/Holliday_junction
    public void ApplyTranslocation(int contigA, int contigB, long posA, long posB, bool inverted)
    {
        var refContig = _contigs[contigA];
        var altContig = _contigs[contigB];
        if (inverted)
        {
            altContig.Invert();
        }
        var splitRef = refContig.Split(posA, true);
        var splitAlt = altContig.Split(posB, true);
        refContig.Join(splitAlt);
        altContig.Join(splitRef);
    }

    public void ApplyWGD()
    {
        _contigs.AddRange(_contigs.Select(ch => new Contig(ch)).ToList());
    }

    public void ApplyChromothripsis(int contigID, List<long> stops, IEnumerable<int> selection)
    {
        var contig = _contigs[contigID];
        contig.ScatterAndGather(stops, selection);
    }
    
    public void ApplyChromoplexy(List<int> contigIDs, List<List<long>> stops, IEnumerable<int> sequence, List<long> breakpoints)
    {
        var subContigs = 
            Enumerable.Range(0, contigIDs.Count)
                .Select(i => _contigs[contigIDs[i]].Scatter(stops[i]))
                .SelectMany(x => x)
                .ToList();
        var ordered = Contig.Concat(sequence.Select(i => subContigs[i]));
        var newContigs = ordered.Scatter(breakpoints).ToList();
        for (int i = 0; i < contigIDs.Count; i++)
        {
            _contigs[contigIDs[i]] = newContigs[i];
        }
    }
    
    public void ApplyPyrgo(int contigID, List<(long start, long len)> frags)
    {
        var contig = _contigs[contigID];
        long offset = 0;
        foreach ((long start, long len) in frags)
        {
            contig.DuplicateRange(start + offset, start + offset + len);
            offset += len;
        }
    }

    // Fragments that do not return to the original chromosome
    public void ApplyTIChain(List<(int id, long start, long len, bool dir)> frags)
    {
        var template = new Contig();
        foreach (var frag in frags)
        {
            template.AppendContig(_contigs[frag.id].GetSubContig(frag.start, frag.start + frag.len));
        }
        _contigs.Add(template);
    }    
    
    // First segment is the host, but there is no repetition
    public void ApplyTIBridge(List<(int id, long start, long len, bool dir)> frags) 
    {
        var host = _contigs[frags[0].id];
        var template = new Contig();
        foreach (var frag in frags.Skip(1))
        {
            var contig = _contigs[frag.id].GetSubContig(frag.start, frag.start + frag.len, frag.dir);
            template.AppendContig(contig);
        }
        host.InsertContig(template, frags[0].start);
    }    
    
    // First segment is the host, with repetition
    public void ApplyTICycle(List<(int id, long start, long len, bool dir)> frags)
    {
        var host = _contigs[frags[0].id];
        var template = new Contig();
        foreach (var frag in frags)
        {
            var contig = _contigs[frag.id].GetSubContig(frag.start, frag.start + frag.len, frag.dir);
            template.AppendContig(contig);
        }
        host.InsertContig(template, frags[0].start);
    }
    
    public void ApplyRigma(int contigID, long rigmaStart, List<long> rigmaLens)
    {
        var contig = _contigs[contigID];
        bool lastWasDeletion = false;
        foreach (long len in rigmaLens)
        {
            if (lastWasDeletion)
            {
                rigmaStart += len;
            }
            else
            {
                contig.DeleteRange(rigmaStart, rigmaStart + len);
            }
            lastWasDeletion = !lastWasDeletion;
        }
    }
    public void ApplyPointMutation(int contigID, long location, Nucleotide newNucleotide)
    {
        var contig = _contigs[contigID];
        contig.PointMutate(location, newNucleotide);
    }
    
    public void MergeRegions()
    {
        _contigs.ForEach(c => c.MergeRegions());
    }
}
