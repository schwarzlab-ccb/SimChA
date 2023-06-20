// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.Computation;
using SimChA.DataTypes;

namespace SimChA.Simulation;

// Note: Empty contigs are retained in the list, but not reported. This way the initial indexing is preserved.
public class Karyotype
{
    public double FitnessVal { get; private set; }
    
    public bool SexXX { get;  }
    
    public int CountContigs() 
        => _contigs.Count(c => c.Any());
    
    public long GenomeLen() 
        => _contigs.Sum(c => c.Length());
    
    public IEnumerable<int> ContigIds() 
        => _contigs.Select((c, i) => (c, i)).Where(t => t.c.Any()).Select(t => t.i);

    private readonly List<Contig> _contigs;
    private readonly Dictionary<ChrNo, List<GenRange>> _missingRanges;
    
    public Karyotype(bool sexXX)
    {
        _contigs = HGRef.GetGenotype(sexXX).Select(region => new Contig(region)).ToList();
        _missingRanges = Enum.GetValues<ChrNo>().ToDictionary(chrNo => chrNo, _ => new List<GenRange>());
        SexXX = sexXX;
    }
    
    public Karyotype(Karyotype other)
    {
        _contigs = other._contigs.Select(ch => new Contig(ch)).ToList();
        _missingRanges = other._missingRanges;
        SexXX = other.SexXX;
    }
    
    public Karyotype(List<Contig> contigs, List<GenRange> missingList, bool sexXX)
    {
        _contigs = contigs;
        _missingRanges = Enum.GetValues<ChrNo>().ToDictionary(chrNo => chrNo, _ => new List<GenRange>());
        foreach (var range in missingList)
        {
            _missingRanges[range.ChrNo].Add(range);
        }
        SexXX = sexXX;
    }

    public override string ToString()
        => CountContigs() > 0 ? "[" + string.Join(";", _contigs.Where(c => c.Any())) + "]" : "[]";

    public void GlueNeighbours()
    {
        foreach (var contig in _contigs)
        {
            contig.GlueNeighbours();
        }
    }
    
    public long MissingLen()
        => _missingRanges.Sum(r => r.Value.Sum(range => range.Length));
    
    public double CalcPloidy()
        => 2.0 * GenomeLen() / HGRef.GetGenomeLen(SexXX);
    
    public double CalcCoverage()
        => (HGRef.GetGenomeLen(SexXX) - MissingLen()) / (double) HGRef.GetGenomeLen(SexXX);
    
    public IEnumerable<Region> FindRegionsOfChr(ChrNo chrNo) 
        => _contigs.SelectMany(c => c.FindRegionsOfChr(chrNo));

    public IList<GenRange> GetMissingOfChr(ChrNo chrNo)
        => _missingRanges[chrNo];

    public static long GetTail(long segLength, Contig contig, bool fiveToThree) 
        => fiveToThree ? segLength : contig.Length() - segLength;
    
    public long ContigLen(int contigId)
        => contigId < _contigs.Count ? _contigs[contigId].Length() : 0;
    
    private static (long start, long end) GetIndices(Contig contig, long position, bool fiveToThree)
        => fiveToThree ? (0, position) : (position, contig.Length());

    public List<Gene> GetPresentGenes(Dictionary<ChrNo, List<Gene>> geneLists)
        => _contigs.SelectMany(c => c.GetPresentGenes(geneLists)).ToList();

    public double UpdateFitness(GenRef geneRef, FitnessParams fParams)
        => FitnessVal = Fitness.Calculate(this, geneRef, fParams);
    
    public void ApplyTailDeletion(int contigID, long tailLen, bool fiveToThree)
    {
        var contig = _contigs[contigID];
        long tailSplit = GetTail(tailLen, contig, fiveToThree);
        (long tailStart, long tailEnd) = GetIndices(contig, tailSplit, fiveToThree);
        contig.DeleteRange(tailStart, tailEnd);
    }

    public void ApplyBFB(int contigID, long tailLen, bool fiveToThree)
    {
        var contig = _contigs[contigID];
        long tailSplit = GetTail(tailLen, contig, fiveToThree);
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
}