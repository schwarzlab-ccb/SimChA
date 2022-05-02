// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

[Serializable]
public struct SimParams
{
    // Simulator
    public int Seed;

    // Experiment
    public uint Repeats;
    public long MaxPop;
    public uint MaxSteps;
    public double CutOff;
    public uint MinPop;

    // Model
    public double Turnover;
    public double Confinement;
    public double MutationProb;
    public double FitnessMean;
    public int InitMut;

    // Function    
    public FitnessAccType FitnessAcc;
    public FitnessSampleType FitnessDist;
    public FitnessEffectType FitnessEffect;
    
    // Output
    public bool Checkpoints;

    // public SimParams()
    // {
    // AberrationRates = Enum.GetValues<AberrationEnum>().ToDictionary(a => a, _ => 1.0);
    // }

    // public Dictionary<AberrationEnum, double> AberrationRates { get; }

    // public double SumRates() => AberrationRates.Sum(ar => ar.Value);
}