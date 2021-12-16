// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.Simulation;

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

    public void ApplyTailDeletion()
    {
        var randomChromosome = RandomChr();
        int chrLength = randomChromosome.Length;
        int deletionLength = _random.Next(0, chrLength);
        (int start, int end) = _random.CoinFlip() ? (0, deletionLength) : (chrLength - deletionLength, chrLength);
        randomChromosome.DeleteRange(start, end);
    }

    public override string ToString()
    {
        return "[\n\n" + string.Join(",\n\t", Chromosomes) + "]\n";
    }
}