// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.DataTypes;

namespace SimChA.Computation;

public static class CellSampling
{
    public static long TotalCount(IEnumerable<SubClone> population)
        => population.Sum(sc => sc.TotalCount);

    public static long AliveCount(IEnumerable<SubClone> population)
        => population.Sum(sc => sc.AliveCount);

    public static long DeadCount(IEnumerable<SubClone> population)
        => population.Sum(sc => sc.DeadCount);
    
    public static long LostCount(IEnumerable<SubClone> population)
        => population.Sum(sc => sc.Lost);

    public static (long, long, long) PopState(List<SubClone> populations) =>
        (TotalCount(populations),  AliveCount(populations), LostCount(populations));
}