// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

[Serializable]
public class SimParams
{
    public Dictionary<AberrationEnum, double> AberrationRates =
        Enum.GetValues<AberrationEnum>().ToDictionary(a => a, a => 1.0);

    public int Seed;
    public bool IsFemale;
    
    public double SumRates()
        => AberrationRates.Sum(ar => ar.Value);
}