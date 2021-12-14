// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public class Karyotype
{
    public List<Chromosome> Chromosomes { get; }

    public Karyotype()
    {
        Chromosomes = ReferenceGenome.Genotype.Select(region => new Chromosome(region)).ToList();
    }
    
    public Karyotype(Karyotype other)
    {
        Chromosomes = new List<Chromosome>(other.Chromosomes);
    }

    public void AddChromosome(Chromosome c)
    {
        Chromosomes.Add(c);
    }
}