// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using MathNet.Numerics.Distributions;

namespace SimChA.DataTypes;

public class Karyotype
{
    private List<Chromosome> Chromosomes { get; set; }

    public int ChromCount => Chromosomes.Count;

    private readonly Random _random = new();

    public Karyotype(bool isFemale)
    {
        var reference = ReferenceGenome.GetGenotype(isFemale).Select(region => new Chromosome(region));
        Chromosomes = new List<Chromosome>(reference);
    }

    public Karyotype(Karyotype other)
    {
        Chromosomes = other.Chromosomes.Select(ch => new Chromosome(ch)).ToList();
    }
    
    // Specifically used in to remove a chromosome lost during division (mis-segregation / BFB)
    public Karyotype(Karyotype other, Chromosome removed)
    {
        Chromosomes = other.Chromosomes
            .Where(ch => ch != removed)
            .Select(ch => new Chromosome(ch))
            .ToList();
    }
    private Chromosome RandomChr()
    {
        return Chromosomes.Shuffle().First();
    }

    private List<Chromosome> RandomChrs(int count)
    {
        return Chromosomes.Shuffle().Take(count).ToList();
    }

    public override string ToString()
    {
        return Chromosomes.Any() ? "[\n\t" + string.Join(",\n\t", Chromosomes) + "\n]\n" : "[]";
    }

    // Get two positions within the chromosome (boundaries are excluded)
    private (int start, int end) GetGammaFraction(Chromosome chr)
    {
        double fraction = Math.Clamp(Gamma.Sample(_random, 1, 1) / 10, 0, 1);
        int segLength = Math.Min((int) (fraction * chr.Length()), chr.Length() - 1);
        int start = DiscreteUniform.Sample(_random, 1, chr.Length() - segLength);
        int end = Math.Min(start + segLength + 1, chr.Length() - 1);
        return (start, end);
    }

    // Get two positions within the chromosome (boundaries are excluded)
    private int GetUniformPosition(Chromosome chr)
    {
        return _random.Next(1, chr.Length() - 1);
    }
    
    public Karyotype ApplyAbberation(AbberationEnum abberation)
    {
        var selectChrs = RandomChrs(2);
        var chr1 = selectChrs[0]; // just a shortcut
        switch (abberation)
        {
            case AbberationEnum.TailDeletion:
                chr1.Split(GetUniformPosition(chr1), _random.CoinFlip());
                return this;

            case AbberationEnum.Missegregation:
                var misKaryotype = new Karyotype(this, chr1);
                Chromosomes.Add(chr1);
                return misKaryotype;

            case AbberationEnum.InternalDuplication:
                (int dupStart, int dupEnd) = GetGammaFraction(chr1);
                chr1.DuplicateRange(dupStart, dupEnd);
                return this;

            case AbberationEnum.InternalDeletion:
                (int delStart, int delEnd) = GetGammaFraction(chr1);
                chr1.DeleteRange(delStart, delEnd);
                return this;

            case AbberationEnum.Translocation:
                var splits
                    = selectChrs.Select(chr
                        => chr.Split(GetUniformPosition(chr), _random.CoinFlip())
                    ).ToList();
                selectChrs[0].Join(splits[1], _random.CoinFlip());
                selectChrs[1].Join(splits[0], _random.CoinFlip());
                return this;

            case AbberationEnum.Inversion:
                (int invStart, int invEnd) = GetGammaFraction(chr1);
                chr1.InvertRange(invStart, invEnd);
                return this;
            
            case AbberationEnum.Duplication:
                var baseKaryotype = new Karyotype(this);
                Chromosomes.Add(chr1);
                return baseKaryotype;

            case AbberationEnum.BreakageFusionBridge:
                var loseKaryotype = new Karyotype(this, chr1);
                chr1.Bridge(GetUniformPosition(chr1), _random.CoinFlip());
                return loseKaryotype;

            case AbberationEnum.Chromothripsis:
                int shardCount = _random.Next(1, (int) Math.Pow(chr1.Length(), 1 / 3f)); // Needs better estimation
                var positions =
                    Enumerable.Range(0, shardCount)
                        .Select(_ => GetUniformPosition(chr1))
                        .Distinct().ToList();
                positions.Sort();
                int count = _random.Next(1, positions.Count);
                chr1.ScatterAndGather(positions, count);
                return this;

            default:
                return this;
        }
    }
}