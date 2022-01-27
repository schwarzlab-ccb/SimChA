// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public struct SimParams
{
    public bool IsFemale;
    public double DivisionRate;
    public double DivisionSlowDown;
    public double MutationRate;
    public double DriverProb;
    public double DeathRate;
    public double SplitRate;
    public double FitnessInc;
    public int InitialPop;

    public readonly Dictionary<AbberationEnum, double> AbberationRates =
        Enum.GetValues<AbberationEnum>().ToDictionary(a => a, a => 1.0);

    public double RatesSum => AbberationRates.Sum(ar => ar.Value);
}