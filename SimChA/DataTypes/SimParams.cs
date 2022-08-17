// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

[Serializable]
public class SimParams
{
    // Simulator
    public int Seed;
    public int StartMut;
    public int StartPop;
    public int SelectionSize;

    // Experiment
    public uint Reps;
    public long CloneTarget;
    public uint MaxSteps;

    // Model
    public double Turnover;
    public double MutationProb;
    public double DeathRate;
    public bool IsFemale;
    
    public Dictionary<AberrationEnum, double> AberrationRates =
        Enum.GetValues<AberrationEnum>().ToDictionary(a => a, a => 1.0);

    public double SumRates()
        => AberrationRates.Sum(ar => ar.Value);
}