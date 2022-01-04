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

    public void Clean()
    {
        Chromosomes.RemoveAll(c => c.Length() <= 0);
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
        return "[\n\t" + string.Join(",\n\t", Chromosomes) + "\n]\n";
    }

    private (int start, int end) GetGammaFraction(Chromosome chr)
    {
        double fraction = Math.Clamp(Gamma.Sample(_random, 1, 1) / 10, 0, 1);
        int segLength = (int) (fraction * chr.Length());
        int start = DiscreteUniform.Sample(_random, 0, chr.Length() - segLength);
        int end = start + segLength + 1;
        return (start, end);
    }

    public Karyotype ApplyAbberation(AbberationEnum abberation)
    {
        var selectChrs = RandomChrs(2);
        var firstChr = selectChrs[0]; // just a shortcut
        switch (abberation)
        {
            case AbberationEnum.TailDeletion:
                int tailLength = _random.Next(0, firstChr.Length());
                firstChr.Split(tailLength, _random.CoinFlip());
                return this;

            case AbberationEnum.Missegregation:
                var misKaryotype = new Karyotype(this);
                Chromosomes.Add(firstChr);
                misKaryotype.Chromosomes.Remove(firstChr);
                return misKaryotype;

            case AbberationEnum.Repetition:
                (int dupStart, int dupEnd) = GetGammaFraction(firstChr);
                firstChr.DuplicateRange(dupStart, dupEnd);
                return this;

            case AbberationEnum.InternalDeletion:
                (int delStart, int delEnd) = GetGammaFraction(firstChr);
                firstChr.DeleteRange(delStart, delEnd);
                return this;

            case AbberationEnum.Translocation:
                var splits
                    = selectChrs.Select(chr
                        => chr.Split(_random.Next(0, chr.Length()), _random.CoinFlip())
                    ).ToList();
                firstChr.Join(splits[1], _random.CoinFlip());
                selectChrs[1].Join(splits[0], _random.CoinFlip());
                return this;

            case AbberationEnum.Inversion:
                (int invStart, int invEnd) = GetGammaFraction(firstChr);
                firstChr.InvertRange(invStart, invEnd);
                return this;

            case AbberationEnum.BreakageFusionBridge:
                var loseKaryotype = new Karyotype(this);
                int breakagePos = _random.Next(0, firstChr.Length());
                firstChr.Bridge(breakagePos, _random.CoinFlip());
                loseKaryotype.Chromosomes.Remove(firstChr);
                return loseKaryotype;

            case AbberationEnum.Chromothripsis:
                int shardCount =
                    _random.Next(1, (int) Math.Pow(firstChr.Length(), 1 / 3f)); // Needs better estimation
                var positions =
                    Enumerable.Range(0, shardCount)
                        .Select(i => _random.Next(1, firstChr.Length() - 1))
                        .Distinct().ToList();
                positions.Sort();
                int count = _random.Next(1, positions.Count);
                firstChr.ScatterAndGather(positions, count);

                return this;

            default:
                return this;
        }
    }
}