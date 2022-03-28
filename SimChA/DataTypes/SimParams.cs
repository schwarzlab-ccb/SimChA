// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

[Serializable]
public struct SimParams
{
    // Simulator
    public int Seed;
    // Experiment
    public uint Repeats;
    public long PopLimit;
    public uint StepLimit;
    public double CutOff;
    public uint InitPop;
    // Model
    public double DivisionRate;
    public double Confinement;
    public double MutationRate;
    public double SplitRate;
    public double FitnessMean;

    // public SimParams()
    // {
    // AberrationRates = Enum.GetValues<AberrationEnum>().ToDictionary(a => a, _ => 1.0);
    // }

    // public Dictionary<AberrationEnum, double> AberrationRates { get; }

    // public double SumRates() => AberrationRates.Sum(ar => ar.Value);
    
}