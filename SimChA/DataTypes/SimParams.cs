// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

[Serializable]
public struct SimParams
{
    public int Seed;
    public int PopLimit;
    public float CutOff;
    public bool IsFemale;
    public double DivisionRate;
    public double DivisionSlowDown;
    public double MutationRate;
    public double DriverProb;
    public double DeathRate;
    public double SplitRate;
    public double FitnessInc;
    public int InitialPop;

    public readonly Dictionary<AberrationEnum, double> AberrationRates =
        Enum.GetValues<AberrationEnum>().ToDictionary(a => a, a => 1.0);

    public double SumRates() => AberrationRates.Sum(ar => ar.Value);
}