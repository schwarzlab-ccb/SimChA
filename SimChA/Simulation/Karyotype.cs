// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using Extreme.Statistics.Distributions;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Misc;

namespace SimChA.Simulation;

// Note: Empty chromosomes are retained in the list, but not reported. This way the initial indexing is preserved.
public class Karyotype
{
    private readonly List<Chromosome> _chromosomes;
    public int ChromCount => _chromosomes.Count(c => c.Any());
    public float FitnessVal { get; private set; }

    public Karyotype(bool isFemale)
    {
        _chromosomes = ReferenceGenome.GetGenotype(isFemale).Select(region => new Chromosome(region)).ToList();
    }

    public Karyotype(Karyotype other)
    {
        _chromosomes = other._chromosomes.Select(ch => new Chromosome(ch)).ToList();
    }

    public override string ToString()
        => ChromCount > 0 ? "[\n\t" + string.Join(",\n\t", _chromosomes.Where(c => c.Any())) + "\n]\n" : "[]";

    public List<Region> GetAllRegions()
    {
        var regions = new List<Region>();
        _chromosomes.ForEach(ch => regions.AddRange(ch.GetAllRegions()));
        return regions;
    }

    // Segment is at most 2 bases shorter than chr
    private static int GetSegLen(Random rnd, Chromosome chr, double meanLen)
    {
        double fraction = Math.Clamp(ExponentialDistribution.Sample(rnd, meanLen), 0, 1);
        return Math.Min((int)Math.Round(fraction * chr.Length()), chr.Length() - 2);
    }

    // Get two positions within the chromosome (boundaries are excluded)
    private static (int start, int end) GetInternalRange(Random rnd, Chromosome chr, double meanLen)
    {
        int segLength = GetSegLen(rnd, chr, meanLen);
        int start = DiscreteUniformDistribution.Sample(rnd, 1, chr.Length() - segLength);
        int end = Math.Min(start + segLength + 1, chr.Length() - 1);
        return (start, end);
    }

    // Get two positions within the chromosome (boundaries are excluded)
    private static int GetUniformPos(Random rnd, Chromosome chr)
        => rnd.Next(1, chr.Length() - 1);

    private static (int, bool) GetTail(Random rnd, Chromosome chr, double meanLen)
    {
        int segLength = GetSegLen(rnd, chr, meanLen);
        bool fromStart = rnd.CoinFlip();
        int pos = fromStart ? segLength - 1 : chr.Length() - segLength - 1;
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

    // Needs better estimation
    private static int GetChromothripsisSiteCount(Random rnd, Chromosome chr)
        => rnd.Next(1, (int)Math.Pow(chr.Length(), 1 / 3f));

    public string ApplyAberration(Random rnd, AberrationEnum aberration, BaseAbbP paramsSet)
    {
        string descriptor;
        var chrIDs = Enumerable.Range(0, _chromosomes.Count).Where(i => _chromosomes[i].Any()).Shuffle(rnd).ToList();
        int chrID = chrIDs[0];
        var chr = _chromosomes[chrID];
        switch (aberration)
        {
            case AberrationEnum.TailDeletion:
                (int tailSplit, bool tailFromStart) = GetTail(rnd, chr, ((FractionAbbP)paramsSet).MeanLength);
                (int tailStart, int tailEnd) = GetIndices(chr, tailSplit, tailFromStart);
                var removed = chr.Split(tailSplit, !tailFromStart);
                descriptor = $"chr:{chrID};start:{tailStart};end{tailEnd}";
                break;

            case AberrationEnum.ChromDeletion:
                descriptor = $"chr:{chrID}";
                chr.Clear();
                break;

            case AberrationEnum.InternalDuplication:
                (int dupStart, int dupEnd) = GetInternalRange(rnd, chr, ((FractionAbbP)paramsSet).MeanLength);
                chr.DuplicateRange(dupStart, dupEnd);
                descriptor = $"chr:{chrID};start:{dupStart};end:{dupEnd}";
                break;

            case AberrationEnum.InternalDeletion:
                (int delStart, int delEnd) = GetInternalRange(rnd, chr, ((FractionAbbP)paramsSet).MeanLength);
                chr.DeleteRange(delStart, delEnd);
                descriptor = $"chr:{chrID};start:{delStart};end:{delEnd}";
                break;

            case AberrationEnum.Translocation:
                // TODO: This is only crossing in the 5-3 direction on both strands. Should be check against literature.
                int altID = chrIDs[1];
                var alt = _chromosomes[altID];
                var splitChr = chr.Split(GetUniformPos(rnd, chr), true);
                var splitAlt = alt.Split(GetUniformPos(rnd, alt), true);
                descriptor = $"chr_A:{chrID};gave:{splitChr.Length()};chr_B:{altID};gave:{splitAlt.Length()};";
                chr.Join(splitAlt);
                alt.Join(splitChr);
                break;

            case AberrationEnum.InternalInversion:
                (int invStart, int invEnd) = GetInternalRange(rnd, chr, ((FractionAbbP)paramsSet).MeanLength);
                chr.InvertRange(invStart, invEnd);
                descriptor = $"chr:{chrID};start:{invStart};end{invEnd}";
                break;

            case AberrationEnum.ChromDuplication:
                _chromosomes.Add(new Chromosome(chr));
                descriptor = $"chr:{chrID}";
                break;

            case AberrationEnum.BreakageFusionBridge:
                (int bfbPos, bool bfbFromStart) = GetTail(rnd, chr, ((FractionAbbP)paramsSet).MeanLength);
                (int bfbStart, int bfbEnd) = GetIndices(chr, bfbPos, bfbFromStart);
                descriptor = $"chr:{chrID};start:{bfbStart};end{bfbEnd}";
                chr.Bridge(bfbPos, bfbFromStart);
                break;

            case AberrationEnum.WholeGenomeDoubling:
                _chromosomes.AddRange(_chromosomes.Select(ch => new Chromosome(ch)).ToList());
                descriptor = "WGD";
                break;

            case AberrationEnum.Chromothripsis:
                int shardCount = GetChromothripsisSiteCount(rnd, chr);
                var stops = Enumerable.Range(0, shardCount).Select(_ => GetUniformPos(rnd, chr)).Distinct().ToList();
                stops.Sort();
                int count = rnd.Next(1, stops.Count);
                int baseCount = chr.Length();
                chr.ScatterAndGather(stops, count, rnd);
                descriptor = $"chr:{chrID};fragments:{stops.Count + 1};lost:{baseCount - chr.Length()}";
                break;

            case AberrationEnum.Chromoplexy:
            default:
                throw new ArgumentOutOfRangeException(nameof(aberration), aberration, null);
        }
        return descriptor;
    }

    private static (int start, int end) GetIndices(Chromosome chr, int position, bool fromStart)
        => fromStart ? (0, position) : (position, chr.Length());

    public List<Gene> GetPresentGenes(Dictionary<ChromNum, List<Gene>> geneLists)
        => _chromosomes.SelectMany(c => c.GetPresentGenes(geneLists)).ToList();
    
    public float UpdateFitness(Dictionary<ChromNum, List<Gene>> essentialGenes, Dictionary<ChromNum, List<Gene>> tsgOgGenes, SimParams simParams)
        => FitnessVal = Fitness.Calculate(this, essentialGenes, tsgOgGenes, simParams);
}