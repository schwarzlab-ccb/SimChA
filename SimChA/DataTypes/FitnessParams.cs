// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

namespace SimChA.DataTypes;

public class FitnessParams
{
    public double Stress { get; }
    public double TsgOg { get; }
    public double Essentiality { get; }
    public double TotalStrength { get; }

    public List<double> ParamsList()
        => new() { Stress, TsgOg, Essentiality, TotalStrength};

    public FitnessParams(double stress, double tsgOg, double essentiality, double totalStrength)
    {
        var sum = stress + tsgOg + essentiality;
        Stress = stress/sum;
        TsgOg = tsgOg/sum;
        Essentiality = 1.0 - Stress - TsgOg;
        TotalStrength = totalStrength;
    }

    public FitnessParams(List<double> fitnessParamList)
    {
        var sum = fitnessParamList[0] + fitnessParamList[1] + fitnessParamList[2];
        Stress = fitnessParamList[0]/sum;
        TsgOg = fitnessParamList[1]/sum;
        Essentiality = 1.0 - Stress - TsgOg;
        TotalStrength = fitnessParamList[3];
    }
}