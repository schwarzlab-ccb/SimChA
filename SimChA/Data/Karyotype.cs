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

    public Dictionary<GeneListType, Dictionary<Gene, int>> GeneCounts { get; private set; }

    public Karyotype(GenRef genRef, SexType sex)
    {
        GenRef = genRef;
        _contigs = genRef.GetGenotype(sex).Select(region => new Contig(region)).ToList();
        Sex = sex;
        GeneCounts = PresentGenes.GetGeneCounts(_contigs);
    }

    public Karyotype(Karyotype other)
    {
        GenRef = other.GenRef;
        _contigs = other._contigs.ConvertAll(ch => new Contig(ch));
        Sex = other.Sex;
        FitnessVal = other.FitnessVal;
        GeneCounts = other.GeneCounts.ToDictionary(
            kvp => kvp.Key,
            kvp => new Dictionary<Gene, int>(kvp.Value)
        );
    }

    public Karyotype(GenRef genRef, List<Contig> contigs, SexType sex)
    {
        GenRef = genRef;
        _contigs = contigs;
        Sex = sex;
        GeneCounts = PresentGenes.GetGeneCounts(_contigs);
    }

    public int CountContigs()
        => _contigs.Count(c => c.Any());

    public long GenomeLen()
        => _contigs.Sum(c => c.Length);

    public IEnumerable<SNV> GetSNVs()
        => _contigs.SelectMany(c => c.GetSNVs()).ToList();

    public Dictionary<string, List<int>> CalcBreaks()
    {
        var breakSets = GenRef.ChrIDsForSex(Sex).ToDictionary(c => c, c => new HashSet<int> {0, GenRef.ChrLengths[c]});
        foreach (var contig in _contigs)
        {
            foreach ((string chrom, var breaks) in contig.CalcBreaks())
            {
                breakSets[chrom].UnionWith(breaks);
            }
        }

        return breakSets.ToDictionary(k => k.Key, v => v.Value.OrderBy(i => i).ToList());
    }

    public List<CopyNumber> CalcCNs(IDictionary<string, List<int>> allBreaks)
        => CopyNumbers.CalcCNs(allBreaks, _contigs);

    public IEnumerable<int> ContigIds()
        => _contigs.Select((c, i) => (c, i)).Where(t => t.c.Any()).Select(t => t.i);

    public override string ToString()
        => CountContigs() > 0 ? "[" + string.Join(";", _contigs.Where(c => c.Any())) + "]" : "[]";

    // @CODY this should be made private, needs to update tests
    public IEnumerable<Region> FindChrRegions(string chrNo)
        => _contigs.SelectMany(c => c.FindChrRegions(chrNo));

    private static long GetTailSplitPos(long segLength, Contig contig, bool fiveToThree)
        => fiveToThree ? segLength : contig.Length - segLength;

    public long ContigLen(int contigId)
        => contigId < _contigs.Count ? _contigs[contigId].Length : 0;

    public List<(long start, long end)> GetCentromeres(int contigId)
        => _contigs[contigId].GetCentromeres(GenRef.Centromeres);

    private static (long start, long end) GetIndices(Contig contig, long position, bool fiveToThree)
        => fiveToThree ? (0, position) : (position, contig.Length);

    // @CODY this should be made private, needs to update test
    public Contig GetContig(int contigID)
        => _contigs[contigID];

    public IEnumerable<string> GetSeq()
        => _contigs.SelectMany(c => c.GetSeq(GenRef).Concat(new List<string> {"\n"}));

    public IEnumerable<string> GetPresentGenes(string chrom, List<Gene> geneList)
        => _contigs.SelectMany(c => c.GetPresentGenes(chrom, geneList));

    public double UpdateFitness(GenRef genRef, FitParams fParams)
        => FitnessVal = Fitness.Calculate(this, genRef, fParams);

    public void ApplyTailDeletion(int contigID, long tailLen, bool fiveToThree)
    {
        var contig = _contigs[contigID];
        var genesBefore = contig.PresentGenes.Genes;
        long tailSplit = GetTailSplitPos(tailLen, contig, fiveToThree);
        (long tailStart, long tailEnd) = GetIndices(contig, tailSplit, fiveToThree);
        contig.DeleteRange(tailStart, tailEnd);
        PresentGenes.UpdateGeneCounts(GeneCounts, genesBefore, contig.PresentGenes.Genes);
    }

    public void ApplyTailDuplication(int contigID, long tailLen, bool fiveToThree)
    {
        var contig = _contigs[contigID];
        long tailSplit = GetTailSplitPos(tailLen, contig, fiveToThree);
        (long tailStart, long tailEnd) = GetIndices(contig, tailSplit, fiveToThree);
        var newTail = new Contig(contig.GetSubContig(tailStart, tailEnd));
        _contigs.Add(newTail);
        PresentGenes.UpdateGeneCounts(GeneCounts, null, newTail.PresentGenes.Genes);
    }

    public void ApplyBFB(int contigID, long tailLen, bool fiveToThree)
    {
        var contig = _contigs[contigID];
        var genesToRemove = contig.PresentGenes.Genes;
        long tailSplit = GetTailSplitPos(tailLen, contig, fiveToThree);
        contig.Bridge(tailSplit, fiveToThree);
        PresentGenes.UpdateGeneCounts(GeneCounts, genesToRemove, contig.PresentGenes.Genes);
    }

    public void ApplyContigDeletion(int contigID)
    {
        var contig = _contigs[contigID];
        var genesToRemove = contig.PresentGenes.Genes;
        contig.Clear();
        PresentGenes.UpdateGeneCounts(GeneCounts, genesToRemove, null);
    }

    public void ApplyContigDuplication(int contigID)
    {
        var contig = _contigs[contigID];
        _contigs.Add(new Contig(contig));
        PresentGenes.UpdateGeneCounts(GeneCounts, null, contig.PresentGenes.Genes);
    }

    public void ApplyInternalDuplication(int contigID, long startPos, long endPos)
    {
        var contig = _contigs[contigID];
        var genesToRemove = contig.PresentGenes.Genes;
        contig.DuplicateRange(startPos, endPos);
        PresentGenes.UpdateGeneCounts(GeneCounts, genesToRemove, contig.PresentGenes.Genes);
    }

    public void ApplyInvertedDuplication(int contigID, long startPos, long endPos)
    {
        var contig = _contigs[contigID];
        var genesToRemove = contig.PresentGenes.Genes;
        contig.DuplicateRange(startPos, endPos);
        contig.InvertRange(endPos, endPos + (endPos - startPos));
        PresentGenes.UpdateGeneCounts(GeneCounts, genesToRemove, contig.PresentGenes.Genes);
    }

    public void ApplyInternalInversion(int contigID, long startPos, long endPos)
    {
        var contig = _contigs[contigID];
        var genesToRemove = contig.PresentGenes.Genes;
        contig.InvertRange(startPos, endPos);
        PresentGenes.UpdateGeneCounts(GeneCounts, genesToRemove, contig.PresentGenes.Genes);
    }

    public void ApplyInternalDeletion(int contigID, long startPos, long endPos)
    {
        var contig = _contigs[contigID];
        var genesToRemove = contig.PresentGenes.Genes;
        contig.DeleteRange(startPos, endPos);
        PresentGenes.UpdateGeneCounts(GeneCounts, genesToRemove, contig.PresentGenes.Genes);
    }

    // TODO: Implement Update GeneCounts
    // Translocation might invert based on the orientation of the holiday Junction https://en.wikipedia.org/wiki/Holliday_junction
    public void ApplyTranslocation(int contigA, int contigB, long posA, long posB, bool inverted)
    {
        var refContig = _contigs[contigA];
        var altContig = _contigs[contigB];
        if (inverted)
        {
            altContig.Revert();
        }

        var splitRef = refContig.Split(posA, true);
        var splitAlt = altContig.Split(posB, true);
        refContig.Join(splitAlt);
        altContig.Join(splitRef);
    }

    public void ApplyWGD()
    {
        _contigs.AddRange(_contigs.Select(ch => new Contig(ch)).ToList());
        PresentGenes.DoubleGeneCounts(GeneCounts);
    }

    public void ApplyChromothripsis(int contigID, List<long> stops, IEnumerable<int> selection)
    {
        var contig = _contigs[contigID];
        var genesToRemove = contig.PresentGenes.Genes;
        contig.ScatterAndGather(stops, selection);
        PresentGenes.UpdateGeneCounts(GeneCounts, genesToRemove, contig.PresentGenes.Genes);
    }

    // TODO: Implement Update GeneCounts
    public void ApplyChromoplexy(List<int> contigIDs, List<List<long>> stops, IEnumerable<int> sequence,
        List<long> breakpoints)
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
        var genesToRemove = contig.PresentGenes.Genes;
        long offset = 0;
        foreach ((long start, long len) in frags)
        {
            contig.DuplicateRange(start + offset, start + offset + len);
            offset += len;
        }
        PresentGenes.UpdateGeneCounts(GeneCounts, genesToRemove, contig.PresentGenes.Genes);
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
        PresentGenes.UpdateGeneCounts(GeneCounts, null, template.PresentGenes.Genes);
    }

    // First segment is the host, but there is no repetition
    public void ApplyTIBridge(List<(int id, long start, long len, bool dir)> frags)
    {
        var host = _contigs[frags[0].id];
        var genesToRemove = host.PresentGenes.Genes;
        var template = new Contig();
        foreach (var frag in frags.Skip(1))
        {
            var contig = _contigs[frag.id].GetSubContig(frag.start, frag.start + frag.len, frag.dir);
            template.AppendContig(contig);
        }

        host.InsertContig(template, frags[0].start);
        PresentGenes.UpdateGeneCounts(GeneCounts, genesToRemove, host.PresentGenes.Genes);
    }

    // First segment is the host, with repetition
    public void ApplyTICycle(List<(int id, long start, long len, bool dir)> frags)
    {
        var host = _contigs[frags[0].id];
        var genesToRemove = host.PresentGenes.Genes;
        var template = new Contig();
        foreach (var frag in frags)
        {
            var contig = _contigs[frag.id].GetSubContig(frag.start, frag.start + frag.len, frag.dir);
            template.AppendContig(contig);
        }

        host.InsertContig(template, frags[0].start);
        PresentGenes.UpdateGeneCounts(GeneCounts, genesToRemove, host.PresentGenes.Genes);
    }

    // TODO: Implement Update GeneCounts
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

    // TODO: Implement Update GeneCounts
    public void ApplyPointMutation(int contigID, long location, Nucleotide newNucleotide)
    {
        var contig = _contigs[contigID];
        var oldNucleotide = GenRef.GetRefBaseFromSeq(contig.GetSeq(GenRef), (int) location);
        contig.PointMutate(location, oldNucleotide, newNucleotide);
    }

    public void MergeRegions()
    {
        _contigs.ForEach(c => c.MergeRegions());
    }
}