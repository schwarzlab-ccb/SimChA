// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.DataTypes;

namespace SimChA.Simulation;

public static class Utility
{
    public static List<uint> CreateCheckpoints(SimParams simParams)
    {
        var checkpoints =new List<uint> { simParams.InitPop };
        while (checkpoints.Last() < simParams.PopLimit)
        {
            checkpoints.Add(checkpoints.Last() * 2);
        }
        return checkpoints;
    }
}