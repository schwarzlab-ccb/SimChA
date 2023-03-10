// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.Misc;

namespace SimChA.Simulation;

// Note: Empty contigs are retained in the list, but not reported. This way the initial indexing is preserved.
public class Karyotype
{
    public double FitnessVal { get; private set; }
    
    public bool SexXX { get;  }

    public string Sex => SexXX ? "XX" : "XY";
    
    public int CountContigs() 
        => _contigs.Count(c => c.Any());
    
    public long CountBases() 
        => _contigs.Sum(c => c.Length());
    
    private readonly List<Contig> _contigs;
    private readonly List<GenRange> _missingRanges;
    
    public Karyotype(bool sexXX)
    {
        _contigs = HGRef.GetGenotype(sexXX).Select(region => new Contig(region)).ToList();
        _missingRanges = new List<GenRange>();
        SexXX = sexXX;
    }
    
    public Karyotype(Karyotype other)
    {
        _contigs = other._contigs.Select(ch => new Contig(ch)).ToList();
        _missingRanges = other._missingRanges;
        SexXX = other.SexXX;
    }
    
    public Karyotype(List<Contig> contigs, List<GenRange> missingRanges, bool sexXX)
    {
        _contigs = contigs;
        _missingRanges = missingRanges;
        SexXX = sexXX;
    }

    public override string ToString()
        => CountContigs() > 0 ? "[" + string.Join(",", _contigs.Where(c => c.Any())) + "]" : "[]";

    public bool IsMissing(GenRange other)
        => _missingRanges.Any(range => range.Overlaps(other));

    public long MissingLen()
        => _missingRanges.Sum(r => r.Length);
    
    public double CalcPloidy()
        => 2.0 * CountBases() / HGRef.GetGenomeLen(SexXX);
    
    public double CalcMissing()
        => 1 - (HGRef.GetGenomeLen(SexXX) - MissingLen()) / (double) HGRef.GetGenomeLen(SexXX);
    
    public IEnumerable<Region> FindRegionsOfChr(ChrNo chrNo) 
        => _contigs.SelectMany(c => c.FindRegionsOfChr(chrNo));

    public static long GetTail(long segLength, Contig contig, bool fiveToThree) 
        => fiveToThree ? segLength : contig.Length() - segLength;
    
    public long ContigLen(int contigId)
        => contigId < _contigs.Count ? _contigs[contigId].Length() : 0;
    
    private static (long start, long end) GetIndices(Contig contig, long position, bool fiveToThree)
        => fiveToThree ? (0, position) : (position, contig.Length());

    public List<Gene> GetPresentGenes(Dictionary<ChrNo, List<Gene>> geneLists)
        => _contigs.SelectMany(c => c.GetPresentGenes(geneLists)).ToList();
    
    public double UpdateFitness(Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> geneLists, FitnessParams fParams)
        => FitnessVal = Fitness.Calculate(this, geneLists, fParams);

    public string ApplyTailDeletion(int contigID, long tailLen, bool fiveToThree)
    {
        var contig = _contigs[contigID];
        long tailSplit = GetTail(tailLen, contig, fiveToThree);
        (long tailStart, long tailEnd) = GetIndices(contig, tailSplit, fiveToThree);
        contig.DeleteRange(tailStart, tailEnd);
        return $"contig:{contigID};start:{tailStart};end{tailEnd}";
    }

    public string ApplyBFB(int contigID, long tailLen, bool fiveToThree)
    {
        var contig = _contigs[contigID];
        long tailSplit = GetTail(tailLen, contig, fiveToThree);
        (long tailStart, long tailEnd) = GetIndices(contig, tailSplit, fiveToThree);
        contig.Bridge(tailSplit, fiveToThree);
        return $"contig:{contigID};start:{tailStart};end{tailEnd}";
    }
    
    public string ApplyContigDeletion(int contigID)
    {
        var contig = _contigs[contigID];
        contig.Clear();
        return $"contig:{contigID}";
    }
    
    public string ApplyContigDuplication(int contigID)
    {
        var contig = _contigs[contigID];
        _contigs.Add(new Contig(contig));
        return $"contig:{contigID}";
    }

    public string ApplyInternalDuplication(int contigID, long startPos, long endPos)
    {
        var contig = _contigs[contigID];
        contig.DuplicateRange(startPos, endPos);
        return  $"contig:{contigID};start:{startPos};end:{endPos}";
    }
    
    public string ApplyInvertedDuplication(int contigID, long startPos, long endPos)
    {
        var contig = _contigs[contigID];
        contig.DuplicateRange(startPos, endPos);
        contig.InvertRange(endPos, endPos + (endPos - startPos));
        return  $"contig:{contigID};start:{startPos};end:{endPos}";
    }
    
    public string ApplyInternalInversion(int contigID, long startPos, long endPos)
    {
        var contig = _contigs[contigID];
        contig.InvertRange(startPos, endPos);
        return  $"contig:{contigID};start:{startPos};end:{endPos}";
    }

    public string ApplyInternalDeletion(int contigID, long startPos, long endPos)
    {
        var contig = _contigs[contigID];
        contig.DeleteRange(startPos, endPos);
        return  $"contig:{contigID};start:{startPos};end:{endPos}";
    }

    // TODO: This is only crossing in the 5-3 direction on both strands. Should be check against literature.
    public string ApplyTranslocation(int contigA, int contigB, long posA, long posB)
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

    public string ApplyChromothripsis(int contigID, List<long> stops, IEnumerable<int> selection)
    {            
        var contig = _contigs[contigID];
        long contigLen = contig.Length();
        contig.ScatterAndGather(stops, selection);
        return $"contig:{contigID};fragments:{stops.Count + 1};lost:{contigLen - contig.Length()}B";
    }
    
    public string ApplyChromoplexy(List<int> contigIDs, List<List<long>> stops, IEnumerable<int> sequence, List<long> breakpoints)
    {
        throw new NotImplementedException();
    }


    public string ApplyAberration(Random rnd, CNEventP cnEventP)
    {
        using var IDsEnumerator = Enumerable
            .Range(0, _contigs.Count)
            .Where(i => _contigs[i].Any())
            .Shuffle(rnd)
            .GetEnumerator();
        IDsEnumerator.MoveNext();
        int contigA = IDsEnumerator.Current;
        long lenA = _contigs[contigA].Length();
        
        switch (cnEventP.Type)
        {
            // Whole chromosome events
            case CNEventType.ChromDeletion:
                return ApplyContigDeletion(contigA);

            case CNEventType.ChromDuplication:
                return ApplyContigDuplication(contigA);
            
            case CNEventType.WholeGenomeDoubling:
                return ApplyWGD();
            
            // Tail events
            case CNEventType.TailDeletion:
            case CNEventType.BreakageFusionBridge:
                long delFraction = Sampling.GetSegLength(rnd, lenA, cnEventP.Params["Mean"]);
                bool delDirection = rnd.CoinFlip();
                return cnEventP.Type == CNEventType.TailDeletion 
                    ? ApplyTailDeletion(contigA, delFraction, delDirection) 
                    : ApplyBFB(contigA, delFraction, delDirection);

            // Internal events
            case CNEventType.InternalDuplication:
            case CNEventType.InternalDeletion:
            case CNEventType.InternalInversion:
            case CNEventType.InvertedDuplication:
                long segLen = Sampling.GetSegLength(rnd, lenA, cnEventP.Params["Mean"]);
                long start = Sampling.GetInternalPos(rnd, lenA- segLen);
                long end = start + segLen;
                return cnEventP.Type switch
                {
                    CNEventType.InternalDuplication => ApplyInternalDuplication(contigA, start, end),
                    CNEventType.InternalDeletion => ApplyInternalDeletion(contigA, start, end),
                    CNEventType.InternalInversion => ApplyInternalInversion(contigA, start, end),
                    CNEventType.InvertedDuplication => ApplyInternalDuplication(contigA, start, end),
                };

            case CNEventType.Translocation:
                IDsEnumerator.MoveNext();
                int contigB = IDsEnumerator.Current;
                long lenB = _contigs[contigB].Length();
                long posA = Sampling.GetInternalPos(rnd, lenA);
                long posB = Sampling.GetInternalPos(rnd, lenB);
                return ApplyTranslocation(contigA, contigB, posA, posB);
            
            case CNEventType.Chromothripsis:
                int shardCount = Sampling.GetChromothripsisSiteCount(rnd, lenA);
                var stops = Sampling.GetStopsForShards(rnd, lenA, shardCount);
                int shardsKept = rnd.Next(1, stops.Count);
                var order = Enumerable.Range(0, shardCount).Shuffle(rnd).Take(shardsKept);
                return ApplyChromothripsis(contigA, stops, order);
   
            case CNEventType.Chromoplexy:
                int chrCount = Sampling.GetChromoplexySiteCount(rnd);
                var contigIDs = new List<int>();
                var stopsForConting = new List<List<long>>();
                long totalLen = 0;
                int totalFrags = 0;
                for (var i = 0; i < chrCount; i++, IDsEnumerator.MoveNext())
                {
                    contigIDs.Add(IDsEnumerator.Current);
                    totalLen += _contigs[IDsEnumerator.Current].Length();
                    int partsCount = Sampling.GetChromothripsisSiteCount(rnd, lenA);
                    totalFrags += partsCount;
                    stopsForConting.Append(Sampling.GetStopsForShards(rnd, lenA, partsCount));
                }
                var sequence = Enumerable.Range(0, totalFrags).Shuffle(rnd);
                var breakpoints = Sampling.GetStopsForShards(rnd, totalLen, chrCount);
                return "";
                
            default:
                throw new ArgumentOutOfRangeException(nameof(cnEventP.Type), cnEventP.Type, null);
        }
    }
}