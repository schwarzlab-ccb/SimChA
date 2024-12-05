// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

namespace SimChA.DataTypes;

public class FitnessParams
{
    public double Stress { get; }
    public double TsgOg { get; }
    public double Essentiality { get; }
    public double TotalStrength { get; }
    public bool Haploinsufficiency { get; }
    private const double EPSILON = 1e-8;

    public bool NormalizeGenes { get; }

    public List<double> ParamsList(bool includeTotalStrength)
        => includeTotalStrength 
           ? new() { Stress, TsgOg, Essentiality, TotalStrength}
           : new() { Stress, TsgOg, Essentiality};

    public FitnessParams(double stress, double tsgOg, double essentiality, double totalStrength, bool haploinsufficiency = false, bool normalizeGenes = false)
    {
        double sum = stress + tsgOg + essentiality;
        if (sum < EPSILON)
        {
          throw new Exception("FitnessParams must have non-zero entries");
        }
	else
	{
	  Stress = stress/sum;
          TsgOg = tsgOg/sum;
          Essentiality = 1.0 - Stress - TsgOg;
	}
        TotalStrength = totalStrength;
        Haploinsufficiency = haploinsufficiency;
        NormalizeGenes = normalizeGenes;
    }
}
