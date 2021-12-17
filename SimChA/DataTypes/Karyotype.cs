// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com
using MathNet.Numerics.Distributions;

namespace SimChA.DataTypes;

public class Karyotype
{
    private LinkedList<Chromosome> Chromosomes { get; }

    private readonly Random _random = new();

    public Karyotype(bool isFemale)
    {
        var reference = ReferenceGenome.GetGenotype(isFemale).Select(region => new Chromosome(region));
        Chromosomes = new LinkedList<Chromosome>(reference);
    }
    
    public Karyotype(Karyotype other)
    {
        Chromosomes = new LinkedList<Chromosome>(other.Chromosomes);
    }

    private Chromosome RandomChr()
    {
        return Chromosomes.Shuffle().First();
    }

    public override string ToString()
    {
        return "[\n\t" + string.Join(",\n\t", Chromosomes) + "]\n";
    }

    
    public Karyotype ApplyAbberation(AbberationEnum abberation)
    {
        var selectChr = RandomChr();
        int chrLength = selectChr.Length();
        switch (abberation)
        {
            case AbberationEnum.TailDeletion:
                int tailLength = _random.Next(0, chrLength);
                (int start, int end) = _random.CoinFlip() ? (0, tailLength) : (chrLength - tailLength, chrLength);
                selectChr.DeleteRange(start, end);
                return this;
            
            case AbberationEnum.Missegregation:
                var pairKaryotype = new Karyotype(this);
                pairKaryotype.Chromosomes.AddLast(selectChr);
                Chromosomes.Remove(selectChr);
                return pairKaryotype;
            
            case AbberationEnum.Duplication:
                Chromosomes.AddLast(new Chromosome(selectChr));
                return this;
                
            case AbberationEnum.InternalDeletion:
                double fraction = Math.Clamp(Gamma.Sample(_random, 1, 1) / 10, 0, 1);
                int delLength = (int)(fraction * chrLength);
                int delStart = DiscreteUniform.Sample(_random, 0, chrLength - delLength);
                selectChr.DeleteRange(delStart, delStart + delLength + 1);
                return this;
                
            case AbberationEnum.Translocation:
            case AbberationEnum.Chromothripsis:
            case AbberationEnum.Inversion:
            case AbberationEnum.BreakageFusionBridge:
            default:
                return this;
        }
    }
}