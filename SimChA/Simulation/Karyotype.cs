// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using Extreme.Statistics.Distributions;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Misc;

namespace SimChA.Simulation;

// Note: Empty contigs are retained in the list, but not reported. This way the initial indexing is preserved.
public class Karyotype
{
    private readonly List<Contig> _contigs;
    public int ContigCount => _contigs.Count(c => c.Any());
    public double FitnessVal { get; private set; }

    public Karyotype(bool isFemale)
    {
        _contigs = ReferenceGenome.GetGenotype(isFemale).Select(region => new Contig(region)).ToList();
    }

    public Karyotype(Karyotype other)
    {
        _contigs = other._contigs.Select(ch => new Contig(ch)).ToList();
    }

    public override string ToString()
        => ContigCount > 0 ? "[\n\t" + string.Join(",\n\t", _contigs.Where(c => c.Any())) + "\n]\n" : "[]";
    
    public IEnumerable<Region> FindRegionsOfChr(ChrNo chrNo) 
        => _contigs.SelectMany(c => c.FindRegionsOfChr(chrNo));

    public static int GetTail(int segLength, Contig contig, bool fiveToThree) 
        => fiveToThree ? segLength - 1 : contig.Length() - segLength - 1;
    
    
    private static (int start, int end) GetIndices(Contig contig, int position, bool fiveToThree)
        => fiveToThree ? (0, position) : (position, contig.Length());

    public List<Gene> GetPresentGenes(Dictionary<ChrNo, List<Gene>> geneLists)
        => _contigs.SelectMany(c => c.GetPresentGenes(geneLists)).ToList();
    
    public double UpdateFitness(Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> geneLists, SimParams simParams)
        => FitnessVal = Fitness.Calculate(this, geneLists, simParams);

    public string ApplyTailDeletion(int contigID, int tailLen, bool fiveToThree)
    {
        var contig = _contigs[contigID];
        int tailSplit = GetTail(tailLen, contig, fiveToThree);
        (int tailStart, int tailEnd) = GetIndices(contig, tailSplit, fiveToThree);
        contig.DeleteRange(tailStart, tailEnd);
        return $"contig:{contigID};start:{tailStart};end{tailEnd}";
    }

    public string ApplyBFB(int contigID, int tailLen, bool fiveToThree)
    {
        var contig = _contigs[contigID];
        int tailSplit = GetTail(tailLen, contig, fiveToThree);
        (int tailStart, int tailEnd) = GetIndices(contig, tailSplit, fiveToThree);
        contig.Bridge(tailSplit, fiveToThree);
        return $"contig:{contigID};start:{tailStart};end{tailEnd}";
    }
    
    public string ApplyChromDeletion(int contigID)
    {
        var contig = _contigs[contigID];
        contig.Clear();
        return $"contig:{contigID}";
    }
    
    public string ApplyChromDuplication(int contigID)
    {
        var contig = _contigs[contigID];
        _contigs.Add(new Contig(contig));
        return $"contig:{contigID}";
    }

    public string ApplyInternalDuplication(int contigID, int startPos, int endPos)
    {
        var contig = _contigs[contigID];
        contig.DuplicateRange(startPos, startPos + endPos);
        return  $"contig:{contigID};start:{startPos};end:{endPos}";
    }
    
    public string ApplyInternalInversion(int contigID, int startPos, int endPos)
    {
        var contig = _contigs[contigID];
        contig.InvertRange(startPos, startPos + endPos);
        return  $"contig:{contigID};start:{startPos};end:{endPos}";
    }

    public string ApplyInternalDeletion(int contigID, int startPos, int endPos)
    {
        var contig = _contigs[contigID];
        contig.DeleteRange(startPos, startPos + endPos);
        return  $"contig:{contigID};start:{startPos};end:{endPos}";
    }

    // TODO: This is only crossing in the 5-3 direction on both strands. Should be check against literature.
    public string ApplyTranslocation(int contigA, int contigB, int posA, int posB)
    {
        var refContig = _contigs[contigA];
        var altContig = _contigs[contigB];
        var splitRef = refContig.Split(posA, true);
        var splitAlt = altContig.Split(posB, true);
        string descriptor = $"contig_A:{contigA};gave:{splitRef.Length()};contig_B:{contigB};gave:{splitAlt.Length()};";
        refContig.Join(splitAlt);
        altContig.Join(splitRef);
        return descriptor;
    }

    public string ApplyWGD()
    {
        _contigs.AddRange(_contigs.Select(ch => new Contig(ch)).ToList());
        return "";
    }

    public string ApplyChromothripsis(int contigID, List<int> stops, IEnumerable<int> selection)
    {            
        var contig = _contigs[contigID];
        int contigLen = contig.Length();
        contig.ScatterAndGather(stops, selection);
        return $"contig:{contigID};fragments:{stops.Count + 1};lost:{contigLen - contig.Length()}B";
    }
    
    public string ApplyAberration(Random rnd, AberrationEnum aberration, BaseAbbP paramsSet)
    {
        using var IDsEnumerator = Enumerable
            .Range(0, _contigs.Count)
            .Where(i => _contigs[i].Any())
            .Shuffle(rnd)
            .GetEnumerator();
        IDsEnumerator.MoveNext();
        int contigA = IDsEnumerator.Current;
        int lenA = _contigs[contigA].Length();
        
        switch (aberration)
        {
            // Whole chromosome events
            case AberrationEnum.ChromDeletion:
                return ApplyChromDeletion(contigA);

            case AberrationEnum.ChromDuplication:
                return ApplyChromDuplication(contigA);
            
            case AberrationEnum.WholeGenomeDoubling:
                return ApplyWGD();
            
            // Tail events
            case AberrationEnum.TailDeletion:
            case AberrationEnum.BreakageFusionBridge:
                int delFraction = Sampling.GetSegLength(rnd, lenA, ((FractionAbbP) paramsSet).MeanLength);
                bool delDirection = rnd.CoinFlip();
                return aberration == AberrationEnum.TailDeletion 
                    ? ApplyTailDeletion(contigA, delFraction, delDirection) 
                    : ApplyBFB(contigA, delFraction, delDirection);

            // Internal events
            case AberrationEnum.InternalDuplication:
            case AberrationEnum.InternalDeletion:
            case AberrationEnum.InternalInversion:
                int segLen = Sampling.GetSegLength(rnd, lenA, ((FractionAbbP) paramsSet).MeanLength);
                int start = Sampling.GetInternalPos(rnd, lenA- segLen);
                int end = start + segLen;
                return aberration switch
                {
                    AberrationEnum.InternalDuplication => ApplyInternalDuplication(contigA, start, end),
                    AberrationEnum.InternalDeletion => ApplyInternalDeletion(contigA, start, end),
                    _ => ApplyInternalInversion(contigA, start, end)
                };

            case AberrationEnum.Translocation:
                IDsEnumerator.MoveNext();
                int contigB = IDsEnumerator.Current;
                int lenB = _contigs[contigB].Length();
                int posA = Sampling.GetInternalPos(rnd, lenA);
                int posB = Sampling.GetInternalPos(rnd, lenB);
                return ApplyTranslocation(contigA, contigB, posA, posB);
            
            case AberrationEnum.Chromothripsis:
                int shardCount = Sampling.GetChromothripsisSiteCount(rnd, lenA);
                var stops = Sampling.GetStopsForShards(rnd, lenA, shardCount);
                int shardsKept = rnd.Next(1, stops.Count);
                var order = Enumerable.Range(0, shardCount).Shuffle(rnd).Take(shardsKept);
                return ApplyChromothripsis(contigA, stops, order);
   
            case AberrationEnum.Chromoplexy:
            default:
                throw new ArgumentOutOfRangeException(nameof(aberration), aberration, null);
        }
    }
}