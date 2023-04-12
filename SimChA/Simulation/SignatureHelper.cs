// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using Extreme.Statistics.Distributions;
using SimChA.DataTypes;

namespace SimChA.Simulation;

public static class SignatureHelper
{
    public static Signature RndSignature(Random rnd, List<Signature> sigs)
    {
        double probSum = sigs.Sum(ev => ev.Prob);
        double sample = ContinuousUniformDistribution.Sample(rnd, 0, probSum);
        foreach (var sig in sigs)
        {
            if (sample <= sig.Prob)
            {
                return sig;
            }
            sample -= sig.Prob;
        }
        // In the case that float-point calculations would cause jumping out of the loop, use the last one
        return sigs.Last();
    }
    
    public static CNEventP RndEventP(Random rnd, List<CNEventP> eventPs)
    {
        double probSum = eventPs.Sum(ev => ev.Prob);
        double sample = ContinuousUniformDistribution.Sample(rnd, 0, probSum);
        foreach (var ev in eventPs)
        {
            if (sample <= ev.Prob)
            {
                return ev;
            }
            sample -= ev.Prob;
        }
        // In the case that float-point calculations would cause jumping out of the loop, use the last one
        return eventPs.Last();
    }
}