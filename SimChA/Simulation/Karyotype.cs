// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using Extreme.Statistics.Distributions;
using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Misc;

namespace SimChA.Simulation;

public class Karyotype
{
    private readonly List<Chromosome> _chromosomes;
    public int ChromCount => _chromosomes.Count;
    public float Fitness { get; private set; }
    
    public Karyotype(bool isFemale)
    {
        _chromosomes = ReferenceGenome.GetGenotype(isFemale).Select(region => new Chromosome(region)).ToList();
    }

    public Karyotype(Karyotype other)
    {
        _chromosomes = other._chromosomes.Select(ch => new Chromosome(ch)).ToList();
    }

    public override string ToString()
        => _chromosomes.Any() ? "[\n\t" + string.Join(",\n\t", _chromosomes) + "\n]\n" : "[]";

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
        return Math.Min((int) Math.Round(fraction * chr.Length()), chr.Length() - 2);
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
        => rnd.Next(1, (int) Math.Pow(chr.Length(), 1 / 3f));

    public string ApplyAberration(Random rnd, AberrationEnum aberration, BaseAbbP paramsSet)
    {
        string descriptor = "";
        int chrID = _chromosomes.GetRandIndex(rnd);
        var chr = _chromosomes[chrID];
        switch (aberration)
        {
            case AberrationEnum.TailDeletion:
                (int delSplit, bool delFromStart) = GetTail(rnd, chr, ((FractionAbbP) paramsSet).MeanLength);
                var removed = chr.Split(delSplit, !delFromStart);
                descriptor = $"chr:{chrID};lost:{removed.Length()}";
                break;

            case AberrationEnum.ChromDeletion:
                descriptor = $"chr:{chrID}";
                chr.Clear();
                break;

            case AberrationEnum.InternalDuplication:
                (int dupStart, int dupEnd) = GetInternalRange(rnd, chr, ((FractionAbbP) paramsSet).MeanLength);
                chr.DuplicateRange(dupStart, dupEnd);
                descriptor = $"chr:{chrID};start:{dupStart};end:{dupEnd}";
                break;

            case AberrationEnum.InternalDeletion:
                (int delStart, int delEnd) = GetInternalRange(rnd, chr, ((FractionAbbP) paramsSet).MeanLength);
                chr.DeleteRange(delStart, delEnd);                
                descriptor = $"chr:{chrID};start:{delStart};end:{delEnd}";
                break;

            case AberrationEnum.Translocation:
                // TODO: This is only crossing in the 5-3 direction on both strands. Should be check against literature.
                var chrIDs = _chromosomes.GetRandIndices( 2, rnd).ToList();
                var chrA = _chromosomes[chrIDs[0]];
                var chrB = _chromosomes[chrIDs[1]];
                var splitA = chrA.Split(GetUniformPos(rnd, chrA), true);
                var splitB = chrB.Split(GetUniformPos(rnd, chrB), true);
                descriptor = $"chr_A:{chrIDs[0]};gave:{splitA.Length()};chr_B:{chrIDs[0]};gave:{splitB.Length()};";
                chrA.Join(splitB);
                chrB.Join(splitA);
                break;

            case AberrationEnum.InternalInversion:
                (int invStart, int invEnd) = GetInternalRange(rnd, chr, ((FractionAbbP) paramsSet).MeanLength);
                chr.InvertRange(invStart, invEnd);
                descriptor = $"chr:{chrID};start:{invStart};end{invEnd}";
                break;

            case AberrationEnum.ChromDuplication:
                _chromosomes.Add(new Chromosome(chr));
                descriptor = "chr:" + chrID;
                break;

            case AberrationEnum.BreakageFusionBridge:
                (int bfbPos, bool bfbFromStart) = GetTail(rnd, chr, ((FractionAbbP) paramsSet).MeanLength);
                chr.Bridge(bfbPos, bfbFromStart);
                descriptor = $"chr:{chrID};position:{bfbPos}{(bfbFromStart ? " from start":"")}";
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
                descriptor =$"chr:{chrID};fragments:{stops.Count+1};lost:{baseCount - chr.Length()}";
                break;

            case AberrationEnum.Chromoplexy:
            default:
                throw new ArgumentOutOfRangeException(nameof(aberration), aberration, null);
        }
        return descriptor;
    }

    private List<Gene> GetPresentGenes(Dictionary<ChromNum, List<Gene>> geneLists)
    {
        List<Gene> presentGenes = new();
        foreach (var chr in _chromosomes)
        {
            foreach (var region in chr.GetAllRegions())
            {
                var chromNum = region.ChromId.ChromNum;
                presentGenes.AddRange(geneLists[chromNum].FindAll(g =>  g.Region.IsInside(region)));
            }
        }
        return presentGenes;
    }
    
    public float UpdateFitness(Dictionary<ChromNum, List<Gene>> essentialGenes, Dictionary<ChromNum, List<Gene>> tsgOgGenes, SimParams simParams)
    {
        float stress = CalcStress(simParams.StressFraction, ChromCount);
        var essentialFound = GetPresentGenes(essentialGenes);
        var tsgOgFound = GetPresentGenes(tsgOgGenes);
        var essentialsMissing = FindMissingGenes(essentialFound, essentialGenes);
        var tsgOgMissing = FindMissingGenes(tsgOgFound, tsgOgGenes);
        var tsgOgCounts = tsgOgFound.GroupBy(x => x).ToList();
        float essentialityFitness = essentialsMissing.Sum(g => g.DeltaFitness);
        // Twice the value for missing genes (ploidy 0), -1 multiplicative factor for each missing gene (ploidy 1),
        // n - 2 for each overrepresented gene (ploidy 2+)
        float tsgOgFitness = 2 * tsgOgMissing.Sum(g => g.DeltaFitness) 
                             -tsgOgCounts.Sum(g => g.Key.DeltaFitness * (g.Count() - 2));
        // parametrized linear combination of factors
        Fitness = stress + simParams.TsgOgFraction * tsgOgFitness + simParams.EssentialFraction * essentialityFitness;
        return Fitness;
    }

    private static List<Gene> FindMissingGenes(List<Gene> presentGenes, Dictionary<ChromNum, List<Gene>> geneList)
    {
        var missingGenes = new List<Gene>();
        foreach (var (_, allGenes) in geneList)
        {
            missingGenes.AddRange(allGenes.Except(presentGenes));
        }
        return missingGenes;
    }

    // Represents the limitation of space in the nucleus - more chromosomes ==> more stress
    // TODO: This needs to be validated
    private static float CalcStress(float stressFactor, int chromCount)
        => stressFactor * (float)Math.Pow(Math.Max(0, chromCount - 46), 2);
}