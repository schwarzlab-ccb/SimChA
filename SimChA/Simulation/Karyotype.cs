// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using Extreme.Statistics.Distributions;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;
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
        => CountContigs() > 0 ? "[" + string.Join(";", _contigs.Where(c => c.Any())) + "]" : "[]";

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

    public string ApplyChainTemplatedInsertions(int contigID, List<Region> regions, int lastContigID)
    {
        var contig = _contigs[contigID];
        regions.AddRange(_contigs[lastContigID].GetRegionsAfterRegion(regions.Last()));
        contig.AddRegions(regions);
        var regionIDs = string.Join(",", regions);
        return $"contig:{contigID};regions:{regionIDs}";
    }    
    
    public string ApplyBridgeTemplatedInsertions(int contigID, List<Region> regions)
    {
        var contig = _contigs[contigID];
        contig.AddRegions(regions);
        var regionIDs = string.Join(",", regions);
        return $"contig:{contigID};regions:{regionIDs}";
    }    
    
    public string ApplyCycleTemplatedInsertions(int contigID, List<Region> regions, Random rnd)
    {
        var contig = _contigs[contigID];
        regions.Add(contig.GetRandomRegion(rnd));
        contig.AddRegions(regions);
        var regionIDs = string.Join(",", regions);
        return $"contig:{contigID};regions:{regionIDs}";
    }


    public string ApplyPyrgo(int contigID, List<(long start, long len)> frags)
    {
        var contig = _contigs[contigID];
        var res = "contig:{contigID};";
        long offset = 0;
        foreach ((long start, long len) in frags)
        {
            contig.DuplicateRange(start + offset, start + offset + len);
            offset += len;
            res += $"start:{start + offset};end:{start + offset + len};";
        }
        return res;
    }

    public static List<(long, long)> GetSubsegments(Random rnd, long start, long fragmentLen, double mean)
    {
        long meanSize = (long) (fragmentLen * mean);
        int fracCount = Sampling.GetSiteCount(rnd, mean);
        var frags = new List<(long, long)>();
        for (int i = 0; i < fracCount; i++)
        {
            long fracLen = Sampling.GetExpSeg(rnd, fragmentLen, meanSize);
            long fracStart = Sampling.GetInternalPos(rnd, fragmentLen - fracLen);
            frags.Add((start + fracStart, fracLen));
        }
        return frags;
    }

    public string ApplyRigma(int contigID, long rigmaStart, List<long> rigmaLens)
    {
        var contig = _contigs[contigID];
        var lastWasDeletion = false;
        var res = "contig:{contigID};";
        foreach (long len in rigmaLens)
        {
            if (lastWasDeletion)
            {
                rigmaStart += len;
            }
            else
            {
                contig.DeleteRange(rigmaStart, rigmaStart + len);
                rigmaStart -= len;
                res += "start:{rigmaStart};end:{rigmaStart + len};";
            }
            lastWasDeletion = !lastWasDeletion;
        }
        return res;
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
                long tailSize = cnEventP.Get("Size", 1_000_000);
                long delFraction = Sampling.GetExpSeg(rnd, lenA, tailSize);
                bool delDirection = rnd.CoinFlip();
                return cnEventP.Type == CNEventType.TailDeletion 
                    ? ApplyTailDeletion(contigA, delFraction, delDirection) 
                    : ApplyBFB(contigA, delFraction, delDirection);

            // Internal events
            case CNEventType.InternalDuplication:
            case CNEventType.InternalDeletion:
            case CNEventType.InternalInversion:
            case CNEventType.InvertedDuplication:
                long internalSize = cnEventP.Get("Size", 100_000);
                long segLen = Sampling.GetExpSeg(rnd, lenA, internalSize);
                long start = Sampling.GetInternalPos(rnd, lenA - segLen);
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
                double invProb = cnEventP.Get("InvProb", 0.0);
                bool inverted = invProb != 0.0 && rnd.CoinFlip(invProb);
                return ApplyTranslocation(contigA, contigB, posA, posB, inverted);
            
            case CNEventType.Chromothripsis:
                long chromothripsisLen = cnEventP.Get("Size", 100_000_000L);
                int shardCount = Sampling.GetSiteCount(rnd, lenA / (double) chromothripsisLen);
                var stops = Sampling.GetStopsForShards(rnd, lenA, shardCount);
                int shardsKept = rnd.Next(1, stops.Count);
                var order = Enumerable.Range(0, shardCount).Shuffle(rnd).Take(shardsKept);
                return ApplyChromothripsis(contigA, stops, order);
            
            case CNEventType.Pyrgo:
                long pyrgoLen = cnEventP.Get("Size", 1_000_000L);
                double pyrgoMean = cnEventP.Get("Mean", 0.1);
                long pyrgoFrag = Sampling.GetExpSeg(rnd, lenA, pyrgoLen);
                long pyrgoStart = Sampling.GetInternalPos(rnd, lenA - pyrgoFrag);
                var frags = GetSubsegments(rnd, pyrgoStart, pyrgoFrag, pyrgoMean);
                return ApplyPyrgo(contigA, frags);
            
            case CNEventType.Rigma:
                long rigmaLen = cnEventP.Get("Size", 1_000_000L);
                double rigmaMean = cnEventP.Get("Mean", 0.1);
                long rigmaStart = Sampling.GetInternalPos(rnd, lenA - rigmaLen);
                int rigmaCount = Sampling.GetSiteCount(rnd, rigmaMean);
                var rigmaStops = Enumerable.Range(0, rigmaCount).Select(i => Sampling.GetExpSeg(rnd, lenA, rigmaMean)).ToList();
                return ApplyRigma(contigA, rigmaStart, rigmaStops);

            case CNEventType.Chromoplexy:
                int chrCount = Sampling.GetChromoplexySiteCount(rnd);
                var contigIDs = new List<int>();
                var stopsForContig = new List<List<long>>();
                var totalLen = 0L;
                var totalFrags = 0;
                for (var i = 0; i < chrCount; i++, IDsEnumerator.MoveNext())
                {
                    contigIDs.Add(IDsEnumerator.Current);
                    long thisLen = _contigs[IDsEnumerator.Current].Length();
                    totalLen += thisLen;
                    int partsCount = Sampling.GetSiteCount(rnd, thisLen / (double) lenA);
                    totalFrags += partsCount;
                    stopsForContig.Add(Sampling.GetStopsForShards(rnd, lenA, partsCount));
                }
                var sequence = Enumerable.Range(0, totalFrags).Shuffle(rnd);
                var breakpoints = Sampling.GetStopsForShards(rnd, totalLen, chrCount);
                return ApplyChromoplexy(contigIDs, stopsForContig, sequence, breakpoints);
            
            case CNEventType.TIChain:
            case CNEventType.TICycle:
            case CNEventType.TIBridge:
                const double probabilityOfSuccess = 0.9; // TODO: should be dependent on the event type possibly,
                                                         // also the value should be better justified
                                                         // see https://www.nature.com/articles/s41586-019-1913-9/figures/9
                int numberOfRegions = GeometricDistribution.Sample(rnd, probabilityOfSuccess);
                var regions = new List<Region>();
                for (var i = 0; i < numberOfRegions; i++, IDsEnumerator.MoveNext())
                {
                    regions.Add(_contigs[IDsEnumerator.Current].GetRandomRegion(rnd));
                }
                return cnEventP.Type switch
                {
                    CNEventType.TIChain => ApplyChainTemplatedInsertions(contigA, regions, IDsEnumerator.Current),
                    CNEventType.TICycle => ApplyCycleTemplatedInsertions(contigA, regions, rnd),
                    CNEventType.TIBridge => ApplyBridgeTemplatedInsertions(contigA, regions)
                };

            default:
                throw new ArgumentOutOfRangeException(nameof(cnEventP.Type), cnEventP.Type, null);
        }
    }


    private string ApplyContigEvent(ContigEventData eventData)
    {
        switch (eventData.EventType)
        {
            // Whole chromosome events
            case CNEventType.ChromDeletion:
                return ApplyContigDeletion(eventData.ContigId);

            case CNEventType.ChromDuplication:
                return ApplyContigDuplication(eventData.ContigId);
            
            default:
                throw new ArgumentOutOfRangeException(nameof(eventData.EventType), eventData.EventType, null);
        }
    }

    private string ApplyTailEvent(TailEventData eventData)
    {
        int contigA = eventData.ContigId;
        long delFraction = eventData.DelFraction;
        bool delDirection = eventData.Direction;
        return eventData.EventType == CNEventType.TailDeletion
            ? ApplyTailDeletion(contigA, delFraction, delDirection)
            : ApplyBFB(contigA, delFraction, delDirection);
    }

    private string ApplyEvent(ContigEventData eventData)
    {
        switch (eventData.EventType)
        {
            // Whole chromosome events
            case CNEventType.ChromDeletion:
                return ApplyContigDeletion(eventData.ContigId);

            case CNEventType.ChromDuplication:
                return ApplyContigDuplication(eventData.ContigId);
            
            default:
                throw new ArgumentOutOfRangeException(nameof(eventData.EventType), eventData.EventType, null);
        }
    }

    private string ApplyEvent(TailEventData eventData)
    {
        int contigA = eventData.ContigId;
        long delFraction = eventData.DelFraction;
        bool delDirection = eventData.Direction;
        return eventData.EventType == CNEventType.TailDeletion
            ? ApplyTailDeletion(contigA, delFraction, delDirection)
            : ApplyBFB(contigA, delFraction, delDirection);
    }
    private string ApplyEvent(ChromothripsisEventData data)
    {
        return ApplyChromothripsis(data.ContigId, data.StopsList, data.GetSelection());
    }
    private string ApplyEvent(ChromoplexyEventData data)
    {
        return ApplyChromoplexy(data.ContigIdList, data.StopsList, data.GetSequence(), data.BreakpointsList);
    }
    private string ApplyEvent(PyrgoEventData data)
    {
        return ApplyPyrgo(data.ContigId, data.FragmentsList);
    }
    private string ApplyEvent(RigmaEventData data)
    {
        return ApplyRigma(data.ContigId, data.Start, data.StopsList);
    }
    private string ApplyEvent(BaseEventData data)
    {
        return ApplyWGD();
    }

    private string ApplyEvent(InternalEventData eventData)
    {
        int contigA = eventData.ContigId;
        long start = eventData.Start;
        long end = eventData.End;
        return eventData.EventType switch
        {
            CNEventType.InternalDuplication => ApplyInternalDuplication(contigA, start, end),
            CNEventType.InternalDeletion => ApplyInternalDeletion(contigA, start, end),
            CNEventType.InternalInversion => ApplyInternalInversion(contigA, start, end),
            CNEventType.InvertedDuplication => ApplyInternalDuplication(contigA, start, end),
        };
    }
    private string ApplyEvent(PairEventData eventData)
    {
        int contigA = eventData.ContigIdList[0];
        int contigB = eventData.ContigIdList[1];
        long posA = eventData.PosA;
        long posB = eventData.PosB;
        bool inverted = eventData.Direction;
        return ApplyTranslocation(contigA, contigB, posA, posB, inverted);
    }
    public string ApplyEventData(dynamic eventData)
    {
        return ApplyEvent(eventData);
        switch (eventData.EventType)
        {
            // Whole chromosome events
            case CNEventType.ChromDeletion:
            case CNEventType.ChromDuplication:
                return ApplyContigEvent(eventData as ContigEventData);
                
            case CNEventType.TailDeletion:
            case CNEventType.BreakageFusionBridge:
                return ApplyTailEvent(eventData as TailEventData);

            case CNEventType.WholeGenomeDoubling:
                return ApplyWGD();

            // Internal events
            case CNEventType.InternalDuplication:
            case CNEventType.InternalDeletion:
            case CNEventType.InternalInversion:
            case CNEventType.InvertedDuplication:
                return ApplyInternalEvent(eventData as InternalEventData);

            case CNEventType.Translocation:
                return ApplyPairEvent(eventData as PairEventData);
            
            case CNEventType.Chromoplexy:
                return ApplyChromoplexyEvent(eventData as ChromoplexyEventData);
            case CNEventType.Chromothripsis:
                return ApplyChromothripsisEvent(eventData as ChromothripsisEventData);
            
            case CNEventType.Pyrgo:
                return ApplyPyrgoEvent(eventData as PyrgoEventData);
            case CNEventType.Rigma:
                return ApplyRigmaEvent(eventData as RigmaEventData);

            default:
                throw new ArgumentOutOfRangeException(nameof(eventData.EventType), eventData.EventType, null);
        }
    }

    private string ApplyPairEvent(PairEventData eventData)
    {
        int contigA = eventData.ContigIdList[0];
        int contigB = eventData.ContigIdList[1];
        long posA = eventData.PosA;
        long posB = eventData.PosB;
        bool inverted = eventData.Direction;
        return ApplyTranslocation(contigA, contigB, posA, posB, inverted);
    }
    
    private string ApplyInternalEvent(InternalEventData eventData)
    {
        int contigA = eventData.ContigId;
        long start = eventData.Start;
        long end = eventData.End;
        return eventData.EventType switch
        {
            CNEventType.InternalDuplication => ApplyInternalDuplication(contigA, start, end),
            CNEventType.InternalDeletion => ApplyInternalDeletion(contigA, start, end),
            CNEventType.InternalInversion => ApplyInternalInversion(contigA, start, end),
            CNEventType.InvertedDuplication => ApplyInternalDuplication(contigA, start, end),
        };
    }

    private string ApplyChromothripsisEvent(ChromothripsisEventData data)
    {
        return ApplyChromothripsis(data.ContigId, data.StopsList, data.GetSelection());
    }

    private string ApplyChromoplexyEvent(ChromoplexyEventData data)
    {
        return ApplyChromoplexy(data.ContigIdList, data.StopsList, data.GetSequence(), data.BreakpointsList);
    }

    private string ApplyPyrgoEvent(PyrgoEventData data)
    {
        return ApplyPyrgo(data.ContigId, data.FragmentsList);
    }

    private string ApplyRigmaEvent(RigmaEventData data)
    {
        return ApplyRigma(data.ContigId, data.Start, data.StopsList);
    }

    public BaseEventData GenerateCNEventProperties(Random rnd, CNEventP cnEventP)
    {
        using var IDsEnumerator = Enumerable
            .Range(0, _contigs.Count)
            .Where(i => _contigs[i].Any())
            .Shuffle(rnd)
            .GetEnumerator();
        IDsEnumerator.MoveNext();
        int contigA = IDsEnumerator.Current;
        long lenA = _contigs[contigA].Length();

        var affectedContigIds = new List<int>();

        switch (cnEventP.Type)
        {
            // Whole chromosome events
            case CNEventType.ChromDeletion:
            case CNEventType.ChromDuplication:
                return new ContigEventData(cnEventP, contigA);
            
            case CNEventType.WholeGenomeDoubling:
                return new BaseEventData(cnEventP);
            
            // Tail events
            case CNEventType.TailDeletion:
            case CNEventType.BreakageFusionBridge:
                long tailSize = cnEventP.Get("Size", 1_000_000);
                long delFraction = Sampling.GetExpSeg(rnd, lenA, tailSize);
                bool delDirection = rnd.CoinFlip();
                return new TailEventData(cnEventP, contigA, delFraction, delDirection);

            // Internal events
            case CNEventType.InternalDuplication:
            case CNEventType.InternalDeletion:
            case CNEventType.InternalInversion:
            case CNEventType.InvertedDuplication:
                long internalSize = cnEventP.Get("Size", 100_000);
                long segLen = Sampling.GetExpSeg(rnd, lenA, internalSize);
                long start = Sampling.GetInternalPos(rnd, lenA - segLen);
                long end = start + segLen;
                return new InternalEventData(cnEventP, contigA, start, end);
            
            case CNEventType.Translocation:
                IDsEnumerator.MoveNext();
                int contigB = IDsEnumerator.Current;
                long lenB = _contigs[contigB].Length();
                long posA = Sampling.GetInternalPos(rnd, lenA);
                long posB = Sampling.GetInternalPos(rnd, lenB);
                double invProb = cnEventP.Get("InvProb", 0.0);
                bool inverted = invProb != 0.0 && rnd.CoinFlip(invProb);
                
                affectedContigIds.Add(contigA);
                affectedContigIds.Add(contigB);
                return new PairEventData(cnEventP, affectedContigIds, posA, posB, inverted);
            
            case CNEventType.Chromothripsis:
                long chromothripsisLen = cnEventP.Get("Size", 100_000_000L);
                int shardCount = Sampling.GetSiteCount(rnd, lenA / (double) chromothripsisLen);
                var stops = Sampling.GetStopsForShards(rnd, lenA, shardCount);
                int shardsKept = rnd.Next(1, stops.Count);
                var order = Enumerable.Range(0, shardCount).Shuffle(rnd).Take(shardsKept).ToList();
                return new ChromothripsisEventData(cnEventP, contigA, stops, order);

            case CNEventType.Chromoplexy:
                int chrCount = Sampling.GetChromoplexySiteCount(rnd);
                var contigIDs = new List<int>();
                var stopsForContig = new List<List<long>>();
                var totalLen = 0L;
                var totalFrags = 0;
                for (var i = 0; i < chrCount; i++, IDsEnumerator.MoveNext())
                {
                    contigIDs.Add(IDsEnumerator.Current);
                    long thisLen = _contigs[IDsEnumerator.Current].Length();
                    totalLen += thisLen;
                    int partsCount = Sampling.GetSiteCount(rnd, thisLen / (double) lenA);
                    totalFrags += partsCount;
                    stopsForContig.Add(Sampling.GetStopsForShards(rnd, lenA, partsCount));
                }
                var sequence = Enumerable.Range(0, totalFrags).Shuffle(rnd).ToList();
                var breakpoints = Sampling.GetStopsForShards(rnd, totalLen, chrCount);
                return new ChromoplexyEventData(cnEventP, contigIDs, stopsForContig, sequence, breakpoints);

            case CNEventType.Pyrgo:
                long pyrgoLen = cnEventP.Get("Size", 1_000_000L);
                double pyrgoMean = cnEventP.Get("Mean", 0.1);
                long pyrgoFrag = Sampling.GetExpSeg(rnd, lenA, pyrgoLen);
                long pyrgoStart = Sampling.GetInternalPos(rnd, lenA - pyrgoFrag);
                var frags = GetSubsegments(rnd, pyrgoStart, pyrgoFrag, pyrgoMean);
                return new PyrgoEventData(cnEventP, contigA, frags);

            case CNEventType.Rigma:
                long rigmaLen = cnEventP.Get("Size", 1_000_000L);
                double rigmaMean = cnEventP.Get("Mean", 0.1);
                long rigmaStart = Sampling.GetInternalPos(rnd, lenA - rigmaLen);
                int rigmaCount = Sampling.GetSiteCount(rnd, rigmaMean);
                var rigmaStops = Enumerable.Range(0, rigmaCount).Select(i => Sampling.GetExpSeg(rnd, lenA, rigmaMean)).ToList();
                return new RigmaEventData(cnEventP, contigA, rigmaStart, rigmaStops);

            default:
                throw new ArgumentOutOfRangeException(nameof(cnEventP.Type), cnEventP.Type, null);
        }
    }
}