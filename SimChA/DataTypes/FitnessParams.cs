// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

namespace SimChA.DataTypes;

public record FitnessParams(double Stress, double TsgOg, double TotalStrength)
{
    public double Essentiality 
        => (Stress + TsgOg <= 1.0) 
            ? 1.0 - (Stress + TsgOg) 
            : throw new Exception("Stress, TsgOg, Essentiality parameters should be Dirichlet randomly distributed variables. Check that Stress and TsgOg sum to <= 1.0.");
};