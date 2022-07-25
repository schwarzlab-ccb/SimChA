// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.DataTypes;

namespace SimChA.Simulation;

public static class Utility
{
    public static List<uint> CreateCheckpoints(uint lowerBound, uint upperBound)
    {

        var checkpoints = new List<uint> { lowerBound };
        while (checkpoints.Last() < upperBound)
        {
            checkpoints.Add(checkpoints.Last() * 2);
        }

        return checkpoints;
    }
}