// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using System.Collections.Generic;
using System.Linq;

namespace SimChA.DataTypes
{
    public class Karyotype
    {
        public List<Chromosome> Chromosomes { get; }

        public Karyotype()
        {
            Chromosomes = ReferenceGenome.Genotype.Select(region => new Chromosome(region)).ToList();
        }

        public void AddChromosome(Chromosome c)
        {
            Chromosomes.Add(c);
        }
    }
}