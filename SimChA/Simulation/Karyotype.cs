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
        string stringIDs = string.Join(",", contigIDs);
        return $"contigs:[{stringIDs}];fragments:{subcontigs.Count}";
    }
    
    public string ApplyPyrgo(int contigID, List<(long start, long len)> frags)
    {
        var contig = _contigs[contigID];
        var res = $"contig:{contigID};";
        long offset = 0;
        foreach ((long start, long len) in frags)
        {
            contig.DuplicateRange(start + offset, start + offset + len);
            offset += len;
            res += $"start:{start + offset};end:{start + offset + len};";
        }
        return res;
    }
    
    private static string DirToStr(bool dir) => dir ? ">" : "<";
    private static string FragsToString(IEnumerable<(int id, long start, long len, bool dir)> frags) 
        => string.Join(",", frags.Select(x => $"({x.id},{x.start},{x.len},{DirToStr(x.dir)})"));

    // Fragments that do not return to the original chromosome
    public string ApplyTIChain(List<(int id, long start, long len, bool dir)> frags)
    {
        var template = new Contig();
        foreach (var frag in frags)
        {
            template.AppendContig(_contigs[frag.id].GetSubContig(frag.start, frag.start + frag.len));
        }
        _contigs.Add(template);
        return FragsToString(frags);
    }    
    
    // First segment is the host, but there is no repetition
    public string ApplyTIBridge(List<(int id, long start, long len, bool dir)> frags) 
    {
        var host = _contigs[frags[0].id];
        var template = new Contig();
        foreach (var frag in frags.Skip(1))
        {
            var contig = _contigs[frag.id].GetSubContig(frag.start, frag.start + frag.len, frag.dir);
            template.AppendContig(contig);
        }
        host.InsertContig(template, frags[0].start);
        return FragsToString(frags);
    }    
    
    // First segment is the host, with repetition
    public string ApplyTICycle(List<(int id, long start, long len, bool dir)> frags)
    {
        var host = _contigs[frags[0].id];
        var template = new Contig();
        foreach (var frag in frags)
        {
            var contig = _contigs[frag.id].GetSubContig(frag.start, frag.start + frag.len, frag.dir);
            template.AppendContig(contig);
        }
        host.InsertContig(template, frags[0].start);
        return FragsToString(frags);
    }
    
    public string ApplyRigma(int contigID, long rigmaStart, List<long> rigmaLens)
    {
        var contig = _contigs[contigID];
        var lastWasDeletion = false;
        var res = $"contig:{contigID};";
        foreach (long len in rigmaLens)
        {
            if (lastWasDeletion)
            {
                rigmaStart += len;
            }
            else
            {
                contig.DeleteRange(rigmaStart, rigmaStart + len);
                res += $"start:{rigmaStart};end:{rigmaStart + len};";
            }
            lastWasDeletion = !lastWasDeletion;
        }
        return res;
    }
    
    public string ApplyCNEvent(Random rnd, CNEventP cnEventP)
    {
        // TODO: Replace for inline in usages
        var eventData = GenerateCNEventData(rnd, cnEventP);
        return eventData.ApplyEvent(this);
    }

    public string ApplyEvent(ContigEventData eventData)
    {
        return eventData.EventType switch
        {
            // Whole chromosome events
            CNEventType.ChromDeletion => ApplyContigDeletion(eventData.ContigId),
            CNEventType.ChromDuplication => ApplyContigDuplication(eventData.ContigId),
            _ => throw new ArgumentOutOfRangeException(nameof(eventData.EventType), eventData.EventType, null)
        };
    }

    public string ApplyEvent(TailEventData eventData)
    {
        int contigA = eventData.ContigId;
        long delFraction = eventData.DelFraction;
        bool delDirection = eventData.Direction;
        return eventData.EventType == CNEventType.TailDeletion
            ? ApplyTailDeletion(contigA, delFraction, delDirection)
            : ApplyBFB(contigA, delFraction, delDirection);
    }
    
    public string ApplyEvent(ChromothripsisEventData data)
    {
        return ApplyChromothripsis(data.ContigId, data.StopsList, data.SelectionList);
    }
    public string ApplyEvent(ChromoplexyEventData data)
    {
        return ApplyChromoplexy(data.ContigIdList, data.StopsList, data.SequenceList, data.BreakpointsList);
    }
    public string ApplyEvent(PyrgoEventData data)
    {
        return ApplyPyrgo(data.ContigId, data.FragmentsList);
    }
    public string ApplyEvent(RigmaEventData data)
    {
        return ApplyRigma(data.ContigId, data.Start, data.StopsList);
    }
    
    public string ApplyEvent(BaseEventData data)
    {
        return ApplyWGD();
    }

    public string ApplyEvent(TemplatedEventData data)
    {
        return data.EventType switch
        {
            CNEventType.TIBridge => ApplyTIBridge(data.Frags),
            CNEventType.TIChain => ApplyTIChain(data.Frags),
            CNEventType.TICycle => ApplyTICycle(data.Frags),
            _ => throw new ArgumentOutOfRangeException(nameof(data.EventType), data.EventType, null)
        };
    }

    public string ApplyEvent(InternalEventData data)
    {
        int contigA = data.ContigId;
        long start = data.Start;
        long end = data.End;
        return data.EventType switch
        {
            CNEventType.InternalDuplication => ApplyInternalDuplication(contigA, start, end),
            CNEventType.InternalDeletion => ApplyInternalDeletion(contigA, start, end),
            CNEventType.InternalInversion => ApplyInternalInversion(contigA, start, end),
            CNEventType.InvertedDuplication => ApplyInternalDuplication(contigA, start, end),
            _ => throw new ArgumentOutOfRangeException(nameof(data.EventType), data.EventType, null)
        };
    }
    public string ApplyEvent(PairEventData eventData)
    {
        long posA = eventData.PosA;
        long posB = eventData.PosB;
        bool inverted = eventData.Direction;
        return ApplyTranslocation(eventData.ContigIdA, eventData.ContigIdB, posA, posB, inverted);
    }

    public BaseEventData GenerateCNEventData(Random rnd, CNEventP cnEventP)
    {
        using var IDsEnumerator = Enumerable
            .Range(0, _contigs.Count)
            .Where(i => _contigs[i].Any())
            .Shuffle(rnd)
            .GetEnumerator();
        IDsEnumerator.MoveNext();
        int contigId = IDsEnumerator.Current;

        switch (cnEventP.Type)
        {
            // Whole chromosome events
            case CNEventType.ChromDeletion:
            case CNEventType.ChromDuplication:
                return new ContigEventData(cnEventP, contigId);
            
            case CNEventType.WholeGenomeDoubling:
                return new BaseEventData(cnEventP);
            
            // Tail events
            case CNEventType.TailDeletion:
            case CNEventType.BreakageFusionBridge:
                return new TailEventData(rnd, this, cnEventP, contigId);

            // Internal events
            case CNEventType.InternalDuplication:
            case CNEventType.InternalDeletion:
            case CNEventType.InternalInversion:
            case CNEventType.InvertedDuplication:
                return new InternalEventData(rnd, this, cnEventP, contigId);
            
            case CNEventType.Translocation:
                IDsEnumerator.MoveNext();
                return new PairEventData(rnd, this, cnEventP, contigId, IDsEnumerator.Current);
            
            case CNEventType.Chromothripsis:
                return new ChromothripsisEventData(rnd, this, cnEventP, contigId);

            case CNEventType.Chromoplexy:
                return new ChromoplexyEventData(rnd, this, cnEventP, IDsEnumerator);

            case CNEventType.Pyrgo:
                return new PyrgoEventData(rnd, this, cnEventP, contigId);

            case CNEventType.Rigma:
                return new RigmaEventData(rnd, this, cnEventP, contigId);
            
            case CNEventType.TIChain:
            case CNEventType.TICycle:
            case CNEventType.TIBridge:
                return new TemplatedEventData(rnd, this, cnEventP, IDsEnumerator);

            default:
                throw new ArgumentOutOfRangeException(nameof(cnEventP.Type), cnEventP.Type, null);
        }
    }
}