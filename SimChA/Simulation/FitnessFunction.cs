// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using MathNet.Numerics.Distributions;
using SimChA.DataTypes;

namespace SimChA.Simulation;

public static class FitnessFunction
{
    public static double SampleFitness(SimParams simParams, Random rnd)
    {
        switch (simParams.FitnessDist)
        {
            case FitnessSampleType.Exponential:
                return Exponential.Sample(rnd, 1/simParams.FitnessMean);
            
            case FitnessSampleType.Beta:
                return Beta.Sample(rnd, 1, 1/simParams.FitnessMean - 1);
                
            case FitnessSampleType.Normal:
                return Math.Max(Normal.Sample(rnd, simParams.FitnessMean, simParams.FitnessMean/2.0), 0);
            
            case FitnessSampleType.Uniform:
                return ContinuousUniform.Sample(rnd, 0, simParams.FitnessMean * 2.0);

            case FitnessSampleType.Constant:
            default:
                return simParams.FitnessMean;
        }
    }
}