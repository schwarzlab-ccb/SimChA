// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using Extreme.Statistics.Distributions;
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
    
    public long GenomeLen() 
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
        => 2.0 * GenomeLen() / HGRef.GetGenomeLen(SexXX);
    
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

    public void GlueNeighbours()
        => _contigs.ForEach(c => c.GlueNeighbours());
    
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

    // Translocation might invert based on the orientation of the holiday Junction https://en.wikipedia.org/wiki/Holliday_junction
    public string ApplyTranslocation(int contigA, int contigB, long posA, long posB, bool inverted)
    {
        var refContig = _contigs[contigA];
        var altContig = _contigs[contigB];
        if (inverted)
        {
            altContig.Invert();
        }
        var splitRef = refContig.Split(posA, true);
        var splitAlt = altContig.Split(posB, true);
        var descriptor = $"contig_A:{contigA};gave:{splitRef.Length()};contig_B:{contigB};gave:{splitAlt.Length()};inverted_B:{inverted}";
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
        var subcontigs = 
            Enumerable.Range(0, contigIDs.Count)
                .Select(i => _contigs[contigIDs[i]].Scatter(stops[i]))
                .SelectMany(x => x)
                .ToList();
        var ordered = Contig.Concat(sequence.Select(i => subcontigs[i]));
        var newContigs = ordered.Scatter(breakpoints).ToList();
        for (var i = 0; i < contigIDs.Count; i++)
        {
            _contigs[contigIDs[i]] = newContigs[i];
        }
        var stringIDs = string.Join(",", contigIDs);
        return $"contigs:[{stringIDs}];fragments:{subcontigs.Count}";
    }

    public string ApplyTIChain(int contigID, List<Region> regions, int lastContigID, long splitPos, Random rnd, double mean)
    {
        var contig = _contigs[contigID];
        regions.AddRange(_contigs[lastContigID].GetSplitRegion(splitPos, regions.Last().Forward));
        splitPos = Sampling.GetSegLength(rnd, contig.Length(), mean);
        contig.AddRegionsChain(regions, splitPos);
        var regionIDs = string.Join(",", regions);
        return $"contig:{contigID};regions:{regionIDs}";
    }    
    
    public string ApplyTIBridge(int contigID, List<Region> regions, Random rnd, double mean)
    {
        var contig = _contigs[contigID];
        var splitPos = Sampling.GetSegLength(rnd, contig.Length(), mean);
        contig.AddRegions(regions, splitPos);
        var regionIDs = string.Join(",", regions);
        return $"contig:{contigID};regions:{regionIDs}";
    }    
    
    public string ApplyTICycle(int contigID, List<Region> regions, Random rnd, double mean)
    {
        var contig = _contigs[contigID];
        var segmentLen = Sampling.GetSegLength(rnd, contig.Length(), mean);
        var startRegion = Sampling.GetInternalPos(rnd, contig.Length());
        var endRegion = startRegion + segmentLen;
        regions.AddRange(contig.GetRandomRegion(startRegion, endRegion));
        contig.Split(startRegion, true);
        var regionIDs = string.Join(",", regions);
        return $"contig:{contigID};regions:{regionIDs}";
    }

    public string ApplyCNEvent(Random rnd, CNEventP cnEventP)
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
                long start = Sampling.GetInternalPos(rnd, lenA - segLen);
                long end = start + segLen;
                return cnEventP.Type switch
                {
                    CNEventType.InternalDuplication => ApplyInternalDuplication(contigA, start, end),
                    CNEventType.InternalDeletion => ApplyInternalDeletion(contigA, start, end),
                    CNEventType.InternalInversion => ApplyInternalInversion(contigA, start, end),
                    CNEventType.InvertedDuplication => ApplyInvertedDuplication(contigA, start, end),
                };

            case CNEventType.Translocation:
                IDsEnumerator.MoveNext();
                int contigB = IDsEnumerator.Current;
                long lenB = _contigs[contigB].Length();
                long posA = Sampling.GetInternalPos(rnd, lenA);
                long posB = Sampling.GetInternalPos(rnd, lenB);
                bool inverted = cnEventP.Params.ContainsKey("InvProb") ? rnd.CoinFlip(cnEventP.Params["InvProb"]) : false;
                return ApplyTranslocation(contigA, contigB, posA, posB, inverted);
            
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
                var totalLen = 0L;
                var totalFrags = 0;
                for (var i = 0; i < chrCount; i++, IDsEnumerator.MoveNext())
                {
                    contigIDs.Add(IDsEnumerator.Current);
                    totalLen += _contigs[IDsEnumerator.Current].Length();
                    int partsCount = Sampling.GetChromothripsisSiteCount(rnd, lenA);
                    totalFrags += partsCount;
                    stopsForConting.Add(Sampling.GetStopsForShards(rnd, lenA, partsCount));
                }
                var sequence = Enumerable.Range(0, totalFrags).Shuffle(rnd);
                var breakpoints = Sampling.GetStopsForShards(rnd, totalLen, chrCount);
                return ApplyChromoplexy(contigIDs, stopsForConting, sequence, breakpoints);
            
            case CNEventType.TIChain:
            case CNEventType.TICycle:
            case CNEventType.TIBridge:
                // TODO: Probability should be dependent on the event type possibly,
                // also the value should be better justified
                // see https://www.nature.com/articles/s41586-019-1913-9/figures/9
                int numberOfRegions = GeometricDistribution.Sample(rnd, cnEventP.Params["Prob"]) + 1;
                var regions = new List<Region>();
                var mean = cnEventP.Params["mean"];
                long endRegion = 0;
                for (var i = 0; i < numberOfRegions; i++, IDsEnumerator.MoveNext())
                {
                    var contigLen = _contigs[IDsEnumerator.Current].Length();
                    var segmentLen = Sampling.GetSegLength(rnd, contigLen, cnEventP.Params["Mean"]);
                    var startRegion = Sampling.GetInternalPos(rnd, contigLen);
                    endRegion = startRegion + segmentLen;
                    regions.AddRange(_contigs[IDsEnumerator.Current].GetRandomRegion(startRegion, endRegion));
                }
                return cnEventP.Type switch
                {
                    CNEventType.TIChain => ApplyTIChain(contigA, regions, IDsEnumerator.Current, endRegion, rnd, mean),
                    CNEventType.TICycle => ApplyTICycle(contigA, regions, rnd, mean),
                    CNEventType.TIBridge => ApplyTIBridge(contigA, regions, rnd, mean)
                };
            
            default:
                throw new ArgumentOutOfRangeException(nameof(cnEventP.Type), cnEventP.Type, null);
        }
    }
}