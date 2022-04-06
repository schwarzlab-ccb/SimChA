// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using MathNet.Numerics.Distributions;
using SimChA.DataTypes;

namespace SimChA.Computation;

public static class CellSampling
{
    public static long PopulationSize(IEnumerable<SubClone> population)
        => population.Sum(sc => sc.TotalCount);

    public static long AliveCount(IEnumerable<SubClone> population)
        => population.Sum(sc => sc.AliveCount);

    
    public static IEnumerable<SubClone> Flatten(IEnumerable<IEnumerable<SubClone>> populations) 
        => populations.SelectMany(x => x);

    public static (long, long) PopState(List<SubClone> populations) =>
        (PopulationSize(populations),  AliveCount(populations));
}