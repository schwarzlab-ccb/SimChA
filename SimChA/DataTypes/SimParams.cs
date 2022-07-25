// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

[Serializable]
public struct SimParams
{
    // Simulator
    public int Seed;
    public int StartMut;
    public int StartPop;

    // Experiment
    public uint Reps;
    public long CloneTarget;
    public uint MaxSteps;

    // Model
    public double Turnover;
    public double MutationProb;
    public double DeathRate;
}