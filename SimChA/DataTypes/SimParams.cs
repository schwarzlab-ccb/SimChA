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
    public double FitnessLambda;
    public double FitnessInc;
    public int InitialPop;
    public bool IsMultiplicative;

    public SimParams()
    {
        Seed = 0;
        PopLimit = 0;
        StepLimit = 0;
        CutOff = 0;
        IsFemale = false;
        DivisionRate = 0;
        Confinement = 0;
        MutationRate = 0;
        DriverProb = 0;
        DeathRate = 0;
        DecayRate = 0;
        SplitRate = 0;
        FitnessLambda = 0;
        FitnessInc = 0;
        InitialPop = 0;
        IsMultiplicative = false;
        AberrationRates = Enum.GetValues<AberrationEnum>().ToDictionary(a => a, _ => 1.0);
    }

    public Dictionary<AberrationEnum, double> AberrationRates { get; }

    public double SumRates() => AberrationRates.Sum(ar => ar.Value);
}