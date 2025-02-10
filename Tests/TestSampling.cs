using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Statistics;
using NUnit.Framework;
using SimChA.Computation;
using SimChA.Simulation;
using SimChA.IO;
using SimChA.DataTypes;

namespace Tests;

[TestFixture]
public class TestSampling
{
    Random _rnd;
    
    [SetUp]
    public void Setup()
    {
        _rnd = new Random(0);
    }

    [Test]
    public void TestDistSampling([Values] DistType dist)
    {
        var reps = Enumerable.Range(0, 100).ToList();
        var res = reps.Select(i => Sampling.SampleDist(_rnd, dist, 100)).Mean();
        Assert.GreaterOrEqual(res, 50);
        Assert.LessOrEqual(res, 150);
    }

    [Test]
    public void TestStops()
    {
        var stops = Sampling.GetStopsForShards(_rnd, 1_000_000, 10);
        Assert.AreEqual(9, stops.Count);
        long prev = 0;
        foreach (long stop in stops)
        {
            Assert.Greater(stop, 0);
            Assert.Less(stop, 1_000_000);
            Assert.Greater(stop, prev);
            prev = stop;
        }
    }
    
    [Test]
    public void TestGetExpSeg()
    {
        Assert.AreEqual(31987, Sampling.GetExpSeg(_rnd, 1_000_000, 0.1));
        Assert.AreEqual(202, Sampling.GetExpSeg(_rnd, 1_000_000, 1_000));
    }
    
    [Test]
    public void TestGetNormSeg()
    {
        Assert.AreEqual(127097, Sampling.GetNormSeg(_rnd, 1_000_000, 0.1));
    }

    [Test]
    public void TestDirichlet()
    {
        var mixture = Sampling.CreateRandomMixture(_rnd, new[] {1000.0, 1.0, 0.001});
        Assert.AreEqual(3, mixture.Count);
        Assert.AreEqual(1.0, mixture.Sum(), 1e-6);
        mixture.ForEach(x => Assert.GreaterOrEqual(x, 0.0));
        Assert.GreaterOrEqual(mixture[0], mixture[1]);
        Assert.GreaterOrEqual(mixture[1], mixture[2]);
    }

    [Test]
    public void TestIndexPicker()
    {
        for (int i = 0; i < 10; i++)
        {
            var rnd = new Random();
            var probs = new List<double> { 0.1, 0.2, 0.0, -1.0, 0.3, 0.4 };
            int index = Extensions.PickRndIndex(rnd, probs);
            Assert.GreaterOrEqual(index, 0);
            Assert.Less(index, probs.Count);
            Assert.AreNotEqual(index, 2);
            Assert.AreNotEqual(index, 3);
        }
    }
}
