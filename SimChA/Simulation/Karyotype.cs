// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.DataTypes;

namespace SimChA.Simulation;

public class Karyotype
{
    private List<Chromosome> Chromosomes { get; }

    private Random Rnd { get; }
    
    public bool IsFemale { get; }

    public int ChromCount => Chromosomes.Count;

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
    
    // Specifically used in to remove a chromosome lost during division (mis-segregation / BFB)
    public Karyotype(Karyotype other, Chromosome removed)
    {
        Chromosomes = other.Chromosomes
            .Where(ch => ch != removed)
            .Select(ch => new Chromosome(ch))
            .ToList();
        IsFemale = other.IsFemale;
        Rnd = other.Rnd;
    }

    public List<Region> GetAllRegions()
    {
        var regions = new List<Region>();
        foreach (var chrom in Chromosomes)
        {
            regions.AddRange(chrom.GetAllRegions());
        }
        return regions;
    }

    private Chromosome RandomChr()
    {
        return Chromosomes.Shuffle(Rnd).First();
    }

    private List<Chromosome> RandomChrs(int count)
    {
        return Chromosomes.Shuffle(Rnd).Take(count).ToList();
    }

    public override string ToString()
    {
        return Chromosomes.Any() ? "[\n\t" + string.Join(",\n\t", Chromosomes) + "\n]\n" : "[]";
    }

    // Get two positions within the chromosome (boundaries are excluded)
    private (int start, int end) GetGammaFraction(Chromosome chr)
    {
        double fraction = Math.Clamp(Extreme.Statistics.Distributions.GammaDistribution.Sample(Rnd, 1, 1) / 10, 0, 1);
        int segLength = Math.Min((int) (fraction * chr.Length()), chr.Length() - 1);
        int start = Extreme.Statistics.Distributions.DiscreteUniformDistribution.Sample(Rnd, 1, chr.Length() - segLength);
        int end = Math.Min(start + segLength + 1, chr.Length() - 1);
        return (start, end);
    }

    // Get two positions within the chromosome (boundaries are excluded)
    private int GetUniformPos(Chromosome chr) 
        => Rnd.Next(1, chr.Length() - 1);

    public Karyotype ApplyAbberation(AberrationEnum aberration)
    {
        var selectChrs = RandomChrs(2);
        var chr1 = selectChrs[0]; // just a shortcut
        switch (aberration)
        {
            case AberrationEnum.TailDeletion:
                chr1.Split(GetUniformPos(chr1), Rnd.CoinFlip());
                return this;

            case AberrationEnum.Missegregation:
                var misKaryotype = new Karyotype(this, chr1);
                Chromosomes.Add(chr1);
                return misKaryotype;

            case AberrationEnum.InternalDuplication:
                (int dupStart, int dupEnd) = GetGammaFraction(chr1);
                chr1.DuplicateRange(dupStart, dupEnd);
                return this;

            case AberrationEnum.InternalDeletion:
                (int delStart, int delEnd) = GetGammaFraction(chr1);
                chr1.DeleteRange(delStart, delEnd);
                return this;

            case AberrationEnum.Translocation:
                var splits
                    = selectChrs.Select(chr
                        => chr.Split(GetUniformPos(chr), Rnd.CoinFlip())
                    ).ToList();
                selectChrs[0].Join(splits[1], Rnd.CoinFlip());
                selectChrs[1].Join(splits[0], Rnd.CoinFlip());
                return this;

            case AberrationEnum.Inversion:
                (int invStart, int invEnd) = GetGammaFraction(chr1);
                chr1.InvertRange(invStart, invEnd);
                return this;
            
            case AberrationEnum.Duplication:
                var baseKaryotype = new Karyotype(this);
                Chromosomes.Add(new Chromosome(chr1));
                return baseKaryotype;

            case AberrationEnum.BreakageFusionBridge:
                var loseKaryotype = new Karyotype(this, chr1);
                chr1.Bridge(GetUniformPos(chr1), Rnd.CoinFlip());
                return loseKaryotype;

            case AberrationEnum.Chromothripsis:
                int shardCount = Rnd.Next(1, (int) Math.Pow(chr1.Length(), 1 / 3f)); // Needs better estimation
                var positions = Enumerable.Range(0, shardCount).Select(_ => GetUniformPos(chr1)).Distinct().ToList();
                positions.Sort();
                int count = Rnd.Next(1, positions.Count);
                chr1.ScatterAndGather(positions, count, Rnd);
                return this;

            case AberrationEnum.WholeGenomeDoubling:
                Chromosomes.AddRange(Chromosomes.Select(ch => new Chromosome(ch)).ToList());
                return this;

            default:
                return this;
        }
    }
}