// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

namespace SimChA.DataTypes;

public class FitnessParams
{
    public double Stress { get; }
    public double TsgOg { get; }
    public double Essentiality { get; }
    public double Delta { get;  }
    public bool Haploinsufficiency { get; }
    
    private const double EPSILON = 1e-8;

    public bool NormalizeGenes { get; }

    public List<double> ParamsList()
        => new() { Stress, TsgOg, Essentiality};

    public FitnessParams(double stress, double tsgOg, double essentiality, double delta, bool haploinsufficiency = false, bool normalizeGenes = false)
    {
        Stress = stress;
	    TsgOg  = tsgOg;
	    Essentiality = essentiality;
        Delta = delta;
        Haploinsufficiency = haploinsufficiency;
        NormalizeGenes = normalizeGenes;
    }
}
