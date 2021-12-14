// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public class Karyotype
{
    private LinkedList<Chromosome> Chromosomes { get; }

    public Karyotype()
    {
        var reference = ReferenceGenome.Genotype.Select(region => new Chromosome(region));
        Chromosomes = new LinkedList<Chromosome>(reference);
    }
    
    public Karyotype(Karyotype other)
    {
        Chromosomes = new LinkedList<Chromosome>(other.Chromosomes);
    }

    public void AddChr(Chromosome chr)
    {
        Chromosomes.AddLast(chr);
    }

    public void RemoveChr(Chromosome chr)
    {
        Chromosomes.Remove(chr);
    }

    public Chromosome RandomChr()
    {
        return Chromosomes.Shuffle().First();
    }

    public Chromosome PopRandomChr()
    {
        var chr = RandomChr();
        Chromosomes.Remove(chr);
        return chr;
    }
}