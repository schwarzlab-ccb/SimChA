// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using Extreme.Statistics.Distributions;
using SimChA.DataTypes;
using SimChA.IO;

namespace SimChA.Simulation;

public class Karyotype
{
    public Karyotype(bool isFemale)
    {
        var reference = ReferenceGenome.GetGenotype(isFemale).Select(region => new Chromosome(region));
        Chromosomes = new List<Chromosome>(reference);
        IsFemale = isFemale;
    }

    public Karyotype(Karyotype other)
    {
        Chromosomes = other.Chromosomes.Select(ch => new Chromosome(ch)).ToList();
        IsFemale = other.IsFemale;
    }

    public Karyotype(List<Chromosome> chromosomes, bool isFemale)
    {
        Chromosomes = chromosomes;
        IsFemale = isFemale;
    }

    public override string ToString()
        => Chromosomes.Any() ? "[\n\t" + string.Join(",\n\t", Chromosomes) + "\n]\n" : "[]";

    private List<Chromosome> Chromosomes { get; }
    public bool IsFemale { get; }
    public int ChromCount => Chromosomes.Count;

    public List<Region> GetAllRegions()
    {
        var regions = new List<Region>();
        Chromosomes.ForEach(ch => regions.AddRange(ch.GetAllRegions()));
        return regions;
    }

    private Chromosome RandomChr(Random rnd)
        => Chromosomes.Shuffle(rnd).First();

    private List<Chromosome> RandomChrs(Random rnd, int count)
        => Chromosomes.Shuffle(rnd).Take(count).ToList();

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

    public void ApplyAberration(Random rnd, AberrationEnum aberration, BaseAbbP paramsSet, string cloneName)
    {
        string region = "";
        var chr = RandomChr(rnd);
        switch (aberration)
        {
            case AberrationEnum.TailDeletion:
                int numberNucleotides = chr._regions[0].End;
                (int delSplit, bool delFromStart) = GetTail(rnd, chr, ((FractionAbbP) paramsSet).MeanLength);
                chr.Split(delSplit, !delFromStart);
                region = $"chromosome:{chr._regions[0].ChromId};nucleotides_deleted:" + 
                    $"{(numberNucleotides - delSplit)}{(delFromStart ? " from start":"")}";
                break;

            case AberrationEnum.ChromDeletion:
                region = $"chromosome:{chr._regions[0].ChromId}";
                Chromosomes.Remove(chr);
                break;

            case AberrationEnum.InternalDuplication:
                (int dupStart, int dupEnd) = GetInternalRange(rnd, chr, ((FractionAbbP) paramsSet).MeanLength);
                chr.DuplicateRange(dupStart, dupEnd);
                region = $"chromosome:{chr._regions[0].ChromId};start:{dupStart};end:{dupEnd}";
                break;

            case AberrationEnum.InternalDeletion:
                (int delStart, int delEnd) = GetInternalRange(rnd, chr, ((FractionAbbP) paramsSet).MeanLength);
                chr.DeleteRange(delStart, delEnd);                
                region = $"chromosome:{chr._regions[0].ChromId};start:{delStart};end:{delEnd}";
                break;

            case AberrationEnum.Translocation:
                // TODO: This is only crossing in the 5-3 direction on both strands. Should be check against literature.
                var chrPair = RandomChrs(rnd, 2);
                var splits = chrPair.Select(c => c.Split(GetUniformPos(rnd, c), true)).ToList();
                chrPair[0].Join(splits[1]);
                chrPair[1].Join(splits[0]);
                region = $"chromosome_1:{chrPair[0]._regions[0].ChromId};added_1:{splits[1]._regions[0].ChromId};start_1:{splits[1]._regions[0].Start};end_1:{splits[1]._regions[0].End}" +
                $"chromosome_2:{chrPair[1]._regions[0].ChromId};added_2:{splits[0]._regions[0].ChromId};start_2:{splits[0]._regions[0].Start};end_2:{splits[0]._regions[0].End}";
                break;

            case AberrationEnum.InternalInversion:
                (int invStart, int invEnd) = GetInternalRange(rnd, chr, ((FractionAbbP) paramsSet).MeanLength);
                chr.InvertRange(invStart, invEnd);
                region = $"chromosome:{chr._regions[0].ChromId};start:{invStart};end{invEnd}";
                break;

            case AberrationEnum.ChromDuplication:
                Chromosomes.Add(new Chromosome(chr));
                region = "chromosome:" + chr._regions[0].ChromId;
                break;

            case AberrationEnum.BreakageFusionBridge:
                (int bfbPos, bool bfbFromStart) = GetTail(rnd, chr, ((FractionAbbP) paramsSet).MeanLength);
                chr.Bridge(bfbPos, bfbFromStart);
                region = $"chromosome:{chr._regions[0].ChromId};position:{bfbPos}{(bfbFromStart ? " from start":"")}";
                break;

            case AberrationEnum.WholeGenomeDoubling:
                Chromosomes.AddRange(Chromosomes.Select(ch => new Chromosome(ch)).ToList());
                break;

            case AberrationEnum.Chromothripsis:
                int shardCount = GetChromothripsisSiteCount(rnd, chr);
                var stops = Enumerable.Range(0, shardCount).Select(_ => GetUniformPos(rnd, chr)).Distinct().ToList();
                stops.Sort();
                int count = rnd.Next(1, stops.Count);
                chr.ScatterAndGather(stops, count, rnd);
                region = $"chromosome:{chr._regions[0].ChromId}";
                foreach (var keptRegion in chr._regions)
                {
                    region += $";start:{keptRegion.Start};end:{keptRegion.End}";
                }
                break;

            default:
                throw new ArgumentOutOfRangeException(nameof(aberration), aberration, null);
        }
        new Abberation(cloneName, aberration.ToString(), region);
    }
}