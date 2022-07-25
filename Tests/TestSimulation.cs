// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using System;
using NUnit.Framework;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.Simulation;
using static Extreme.Statistics.Distributions.BinomialDistribution;

namespace Tests;

[TestFixture]
public class TestSimulation
{
    private SimParams _simParams;
    private Simulator _sim;
    
    [SetUp]
    public void Setup()
    {
        _simParams = new SimParams
        {
            Turnover = 0.1f, 
            MutationProb = 0.01f
        };
        _sim = new Simulator(_simParams, new Random(0));
    }

    [Test]
    public void TestBinomial()
    {
        var rnd1 = new Random(0);
        var rnd2 = new Random(0);
        int res1 = 0, res2 = 0;
        for (int i = 0; i < 100; i++)
        {
            res1 = Sample(rnd1, 1000, .5f);
        }
        for (int i = 0; i < 100; i++)
        {
            res2 = Sample(rnd2, 1000, .5f);
        }
        Assert.AreEqual(res1, res2);
    }
}