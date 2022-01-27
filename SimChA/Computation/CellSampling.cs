// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using MathNet.Numerics.Distributions;
using SimChA.DataTypes;

namespace SimChA.Computation;

public static class CellSampling
{
    public static List<SubClone> SampleCells(List<SubClone> population, long sampleSize)
    {
        throw new NotImplementedException();

        // If re-implemented this should be seeded
        Random rnd = new();
        List<SubClone> sample = new();
            
        long popSize = PopulationSize(population);
        double distRatio = Math.Clamp((double) sampleSize / popSize, 0.0, 1.0);

        foreach (var subClone in population)
        {
            int sampled = Binomial.Sample(rnd, distRatio, subClone.AliveCount);
            if (sampled > 0)
            {
                sample.Add(new SubClone(subClone)
                {
                     
                });
            }
        }

        return sample;
    }

    public static long PopulationSize(IEnumerable<SubClone> population)
        => population.Sum(sc => sc.AliveCount);
    
    public static long PopulationSize(IEnumerable<IEnumerable<SubClone>> populations)
        => populations.Sum(PopulationSize);

    public static IEnumerable<SubClone> Flatten(IEnumerable<IEnumerable<SubClone>> populations) 
        => populations.SelectMany(x => x);
}