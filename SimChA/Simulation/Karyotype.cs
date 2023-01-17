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

    // Segment is at most 2 bases shorter than contig
    private static int GetSegLength(Random rnd, Contig contig, double meanLen)
    {
        double fraction = Math.Clamp(ExponentialDistribution.Sample(rnd, meanLen), 0, 1);
        return Math.Min((int)Math.Round(fraction * contig.Length()), contig.Length() - 2);
    }

    // Get two positions within the contig (boundaries are excluded)
    private static (int start, int end) GetInternalRange(Random rnd, Contig contig, double meanLen)
    {
        int segLength = GetSegLength(rnd, contig, meanLen);
        int start = DiscreteUniformDistribution.Sample(rnd, 1, contig.Length() - segLength);
        int end = Math.Min(start + segLength + 1, contig.Length() - 1);
        return (start, end);
    }

    // Get two positions within the contig (boundaries are excluded)
    private static int GetUniformPos(Random rnd, Contig contig)
        => rnd.Next(1, contig.Length() - 1);

    private static (int, bool) GetTail(Random rnd, Contig contig, double meanLen)
    {
        int segLength = GetSegLength(rnd, contig, meanLen);
        bool fromStart = rnd.CoinFlip();
        int pos = fromStart ? segLength - 1 : contig.Length() - segLength - 1;
        return (pos, fromStart);
    }

    // https://ashpublications.org/blood/article/134/Supplement_1/3767/424006/Chromoplexy-and-Chromothripsis-Are-Important
    private static int GetChromoplexySiteCount(Random rnd)
        => rnd.NextSingle() switch
        {
            var n when n < .46 => 3,
            var n when n < .64 => 4,
            var n when n < .74 => 5,
            var n when n < .79 => 6,
            _ => 2
        };

    
    public string ApplyAberration(Random rnd, AberrationEnum aberration, BaseAbbP paramsSet)
    {
        string descriptor;
        var contigIDs = Enumerable.Range(0, _contigs.Count).Where(i => _contigs[i].Any()).Shuffle(rnd).ToList();
        int contigID = contigIDs[0];
        var contig = _contigs[contigID];
        switch (aberration)
        {
            case AberrationEnum.TailDeletion:
                (int tailSplit, bool fiveToThree) = GetTail(rnd, contig, ((FractionAbbP)paramsSet).MeanLength);
                (int tailStart, int tailEnd) = GetIndices(contig, tailSplit, fiveToThree);
                descriptor = $"contig:{contigID};start:{tailStart};end{tailEnd}";
                break;

            case AberrationEnum.ChromDeletion:
                descriptor = $"contig:{contigID}";
                contig.Clear();
                break;

            case AberrationEnum.InternalDuplication:
                (int dupStart, int dupEnd) = GetInternalRange(rnd, contig, ((FractionAbbP)paramsSet).MeanLength);
                contig.DuplicateRange(dupStart, dupEnd);
                descriptor = $"contig:{contigID};start:{dupStart};end:{dupEnd}";
                break;

            case AberrationEnum.InternalDeletion:
                (int delStart, int delEnd) = GetInternalRange(rnd, contig, ((FractionAbbP)paramsSet).MeanLength);
                contig.DeleteRange(delStart, delEnd);
                descriptor = $"contig:{contigID};start:{delStart};end:{delEnd}";
                break;

            case AberrationEnum.Translocation:
                // TODO: This is only crossing in the 5-3 direction on both strands. Should be check against literature.
                int altID = contigIDs[1];
                var alt = _contigs[altID];
                var splitContig = contig.Split(GetUniformPos(rnd, contig), true);
                var splitAlt = alt.Split(GetUniformPos(rnd, alt), true);
                descriptor = $"contig_A:{contigID};gave:{splitContig.Length()};contig_B:{altID};gave:{splitAlt.Length()};";
                contig.Join(splitAlt);
                alt.Join(splitContig);
                break;

            case AberrationEnum.InternalInversion:
                (int invStart, int invEnd) = GetInternalRange(rnd, contig, ((FractionAbbP)paramsSet).MeanLength);
                contig.InvertRange(invStart, invEnd);
                descriptor = $"contig:{contigID};start:{invStart};end{invEnd}";
                break;

            case AberrationEnum.ChromDuplication:
                _contigs.Add(new Contig(contig));
                descriptor = $"contig:{contigID}";
                break;

            case AberrationEnum.BreakageFusionBridge:
                (int bfbPos, bool bfbFromStart) = GetTail(rnd, contig, ((FractionAbbP)paramsSet).MeanLength);
                (int bfbStart, int bfbEnd) = GetIndices(contig, bfbPos, bfbFromStart);
                descriptor = $"contig:{contigID};start:{bfbStart};end{bfbEnd}";
                contig.Bridge(bfbPos, bfbFromStart);
                break;

            case AberrationEnum.WholeGenomeDoubling:
                _contigs.AddRange(_contigs.Select(ch => new Contig(ch)).ToList());
                descriptor = "WGD";
                break;

            case AberrationEnum.Chromothripsis:
                int shardCount = Sampling.GetChromothripsisSiteCount(rnd, contig);
                var stops = Enumerable.Range(0, shardCount).Select(_ => GetUniformPos(rnd, contig)).Distinct().ToList();
                stops.Sort();
                int count = rnd.Next(1, stops.Count);
                int baseCount = contig.Length();
                contig.ScatterAndGather(stops, count, rnd);
                descriptor = $"contig:{contigID};fragments:{stops.Count + 1};lost:{baseCount - contig.Length()}";
                break;

            case AberrationEnum.Chromoplexy:
            default:
                throw new ArgumentOutOfRangeException(nameof(aberration), aberration, null);
        }
        return descriptor;
    }

    private static (int start, int end) GetIndices(Contig contig, int position, bool fiveToThree)
        => fiveToThree ? (0, position) : (position, contig.Length());

    public List<Gene> GetPresentGenes(Dictionary<ChrNo, List<Gene>> geneLists)
        => _contigs.SelectMany(c => c.GetPresentGenes(geneLists)).ToList();
    
    public double UpdateFitness(
        Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> geneLists, 
        SimParams simParams)
        => FitnessVal = Fitness.Calculate(this, geneLists, simParams);
}