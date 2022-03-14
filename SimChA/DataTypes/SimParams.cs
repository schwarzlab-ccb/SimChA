// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

[Serializable]
public struct SimParams
{
    public int Seed;
    public int PopLimit;
    public int StepLimit;
    public float CutOff;
    public bool IsFemale;
    public double DivisionRate;
    public double Confinement;
    public double MutationRate;
    public double DriverProb;
    public double DeathRate;
    public double DecayRate;
    public double SplitRate;
    public double FitnessIncMu;
    public double FitnessIncSigma;
    public int InitialPop;
    public bool IsMultiplicative;

    public readonly Dictionary<AberrationEnum, double> AberrationRates =
        Enum.GetValues<AberrationEnum>().ToDictionary(a => a, a => 1.0);

    public double SumRates() => AberrationRates.Sum(ar => ar.Value);
}