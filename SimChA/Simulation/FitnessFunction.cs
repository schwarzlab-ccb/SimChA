// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using MathNet.Numerics.Distributions;
using SimChA.DataTypes;

namespace SimChA.Simulation;

public static class FitnessFunction
{
    public static double SampleFitness(SimParams simParams, ProgramConfig config, Random rnd)
    {
        switch (config.FitnessType)
        {
            
            case FitnessSampleType.Exponential:
                return Math.Min(Exponential.Sample(rnd, 1/simParams.FitnessMean) * simParams.DivisionRate, 3 * simParams.FitnessMean);
            
            // TODO 
            case FitnessSampleType.Beta:
                return Beta.Sample(rnd, 1, simParams.FitnessMean) * simParams.DivisionRate;
                
            // TODO 
            case FitnessSampleType.Normal:
                return Math.Max(Normal.Sample(simParams.FitnessMean, simParams.FitnessMean) * simParams.DivisionRate, 0);
                
                
            case FitnessSampleType.Constant:
            default:
                return simParams.FitnessMean * simParams.DivisionRate;
        }
    }
}