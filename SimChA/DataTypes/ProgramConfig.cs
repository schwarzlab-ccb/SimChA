// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using System.Reflection.Metadata;

namespace SimChA.DataTypes;

public enum FitnessSampleType { Constant, Normal, Exponential, Beta, Uniform }

public struct ProgramConfig
{
    public bool MultiplicativeFitness;
    public bool StochasticCellLife;
    public FitnessSampleType FitnessType;
}