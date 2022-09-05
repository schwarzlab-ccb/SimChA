// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using Extreme.Statistics.Distributions;
using SimChA.DataTypes;
using SimChA.IO;

namespace SimChA.Simulation;

public class Karyotype
{
    public Karyotype(bool isFemale, Random random)
    {
        var reference = ReferenceGenome.GetGenotype(isFemale).Select(region => new Chromosome(region));
        Chromosomes = new List<Chromosome>(reference);
        IsFemale = isFemale;
        Rnd = random;
    }

    public Karyotype(Karyotype other)
    {
        Chromosomes = other.Chromosomes.Select(ch => new Chromosome(ch)).ToList();
        IsFemale = other.IsFemale;
        Rnd = other.Rnd;
    }
    
    public Karyotype(List<Chromosome> chromosomes, Random rnd, bool isFemale)
    {
        Chromosomes = chromosomes;
        Rnd = rnd;
        IsFemale = isFemale;
    }
    
    public override string ToString()
        => Chromosomes.Any() ? "[\n\t" + string.Join(",\n\t", Chromosomes) + "\n]\n" : "[]";
    
    private List<Chromosome> Chromosomes { get; }
    private Random Rnd { get; }
    public bool IsFemale { get; }
    public int ChromCount => Chromosomes.Count;

    public List<Region> GetAllRegions()
    {
        var regions = new List<Region>();
        Chromosomes.ForEach(ch => regions.AddRange(ch.GetAllRegions()));
        return regions;
    }

    private Chromosome RandomChr() => Chromosomes.Shuffle(Rnd).First();

    private List<Chromosome> RandomChrs(int count) => Chromosomes.Shuffle(Rnd).Take(count).ToList();

    // Segment is at most 2 bases shorter than chr
    private int GetSegLen(Chromosome chr, double meanLen)
    {
        double fraction = ExponentialDistribution.Sample(Rnd, 1 / meanLen);
        return Math.Min((int)(fraction * chr.Length()), chr.Length() - 2); 
    }
    
    // Get two positions within the chromosome (boundaries are excluded)
    private (int start, int end) GetInternalRange(Chromosome chr, double meanLen)
    {
        int segLength = GetSegLen(chr, meanLen);   
        int start = DiscreteUniformDistribution.Sample(Rnd, 1, chr.Length() - segLength);
        int end = Math.Min(start + segLength + 1, chr.Length() - 1);
        return (start, end);
    }

    // Get two positions within the chromosome (boundaries are excluded)
    private int GetUniformPos(Chromosome chr)
        => Rnd.Next(1, chr.Length() - 1);
    
    private (int, bool) GetTail(Chromosome chr, double meanLen)
    {
        int segLength = GetSegLen(chr, meanLen);
        bool fromStart = Rnd.CoinFlip();
        int pos = fromStart ? segLength - 1 : chr.Length() - segLength - 1;
        return (pos, fromStart);
    }

    // https://ashpublications.org/blood/article/134/Supplement_1/3767/424006/Chromoplexy-and-Chromothripsis-Are-Important
    private int GetChromoplexySiteCount()
        => Rnd.NextSingle() switch
        {
            var n when n < .46 => 3,
            var n when n < .64 => 4,
            var n when n < .74 => 5,
            var n when n < .79 => 6,
            _ => 2
        };

    // Needs better estimation
    private int GetChromothripsisSiteCount(Chromosome chr) 
        => Rnd.Next(1, (int)Math.Pow(chr.Length(), 1 / 3f)); 

    public void ApplyAberration(AberrationEnum aberration, BaseAbbP p_abber)
    {
        var chr = RandomChr();
        switch (aberration)
        {
            case AberrationEnum.TailDeletion:
                (int delSplit, bool delFromStart) = GetTail(chr, p_abber.Likelihood);
                chr.Split(delSplit, delFromStart);
                break;

            case AberrationEnum.ChromDeletion:
                Chromosomes.Remove(chr);
                break;
            
            case AberrationEnum.InternalDuplication:
                (int dupStart, int dupEnd) = GetInternalRange(chr, p_abber.Likelihood);
                chr.DuplicateRange(dupStart, dupEnd);
                break;

            case AberrationEnum.InternalDeletion:
                (int delStart, int delEnd) = GetInternalRange(chr, p_abber.Likelihood);
                chr.DeleteRange(delStart, delEnd);
                break;

            case AberrationEnum.Translocation:
                var chrPair = RandomChrs(2);
                var splits = chrPair
                    .Select(c => c.Split(GetUniformPos(c), Rnd.CoinFlip()))
                    .ToList();
                chrPair[0].Join(splits[1], Rnd.CoinFlip());
                chrPair[1].Join(splits[0], Rnd.CoinFlip());
                break;

            case AberrationEnum.Inversion:
                (int invStart, int invEnd) = GetInternalRange(chr, p_abber.Likelihood);
                chr.InvertRange(invStart, invEnd);
                break;

            case AberrationEnum.ChromDuplication:
                Chromosomes.Add(new Chromosome(chr));
                break;

            case AberrationEnum.BreakageFusionBridge:
                (int bfbPos, bool bfbFromStart) = GetTail(chr, p_abber.Likelihood);
                chr.Bridge(bfbPos, bfbFromStart);
                break;

            case AberrationEnum.WholeGenomeDoubling:
                Chromosomes.AddRange(Chromosomes.Select(ch => new Chromosome(ch)).ToList());
                break;
            
            case AberrationEnum.Chromothripsis:
                int shardCount = GetChromothripsisSiteCount(chr);
                var stops = Enumerable.Range(0, shardCount).Select(_ => GetUniformPos(chr)).Distinct().ToList();
                stops.Sort();
                int count = Rnd.Next(1, stops.Count);
                chr.ScatterAndGather(stops, count, Rnd);
                break;
            
            default:
                throw new ArgumentOutOfRangeException(nameof(aberration), aberration, null);
        }
    }
}