// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.DataTypes;

namespace SimChA.Computation;

public static class CellSampling
{
    public static long TotalCount(IEnumerable<SubClone> population)
        => population.Sum(sc => sc.TotalCount);

    public static long AliveCount(IEnumerable<SubClone> population)
        => population.Sum(sc => sc.AliveCount);

    public static long NecroCount(IEnumerable<SubClone> population)
        => population.Sum(sc => sc.NecroCount);
    
    public static long LostCount(IEnumerable<SubClone> population)
        => population.Sum(sc => sc.LostCount);

    public static (long, long, long) PopState(List<SubClone> populations) =>
        (TotalCount(populations),  AliveCount(populations), NecroCount(populations));
}