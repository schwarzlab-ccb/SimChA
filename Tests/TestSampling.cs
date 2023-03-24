using System;
using System.Collections.Generic;
using NUnit.Framework;
using SimChA.Simulation;

namespace Tests;

[TestFixture]

public class TestSampling
{

    [Test]
    public void TestBetaDistribution()
    {
        var rnd = new Random();
        var distribution = new Dictionary<int, int>();
        for(int i = 2; i <= 7; i++)
        {
            distribution.Add(i, 0);
        }
        for(int i = 0; i < 100000; i++)
        {
            distribution[Sampling.BetaDistribution(rnd)]++;
        }

        Assert.GreaterOrEqual(distribution[2], distribution[3]);
        Assert.GreaterOrEqual(distribution[3], distribution[4]);
        Assert.GreaterOrEqual(distribution[4], distribution[5]);
        Assert.GreaterOrEqual(distribution[5], distribution[6]);
        Assert.GreaterOrEqual(distribution[6], distribution[7]);
        Console.WriteLine(distribution);
    
    }
}