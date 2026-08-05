using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Statistics;
using NUnit.Framework;
using SimChA.Computation;
using SimChA.Simulation;

namespace Tests;

[TestFixture]
public class TestSampling
{
    private Random _rnd = null!;
    
    [SetUp]
    public void Setup()
    {
        _rnd = new Random(0);
    }

    [TestCase(0.001), TestCase(0.01), TestCase(0.1), TestCase(0.5)]
    public void TestParetoSampling(double mean)
    {
        double res = Enumerable.Range(0, 1000).Select(i => Sampling.SampleParetoLim(_rnd, mean)).Mean();
        Assert.Greater(res, mean * 0.5);
        Assert.Less(res, mean * 1.5);
        Console.Write(res);
    }

    [Test]
    public void TestContSampling([Values(DistType.Exponential, DistType.Normal)] DistType dist, [Values(0.01, 0.1, 1, 10, 100)] double mean)
    {
        double res = Enumerable.Range(0, 1000).Select(i => Sampling.SampleContDist(_rnd, dist, mean)).Average();
        Assert.Greater(res, mean * 0.5);
        Assert.Less(res, mean * 1.5);
    }
    
    [Test]
    public void TestDiscSampling([Values(DistType.Geometric, DistType.Poisson)] DistType dist, [Values(1, 10, 100)] double mean)
    {
        double res = Enumerable.Range(0, 1000).Select(i => Sampling.SampleDiscDist(_rnd, dist, mean)).Average();
        Assert.Greater(res, mean * 0.5);
        Assert.Less(res, mean * 1.5);
    }

    [Test]
    public void TestGeometric()
    {
        double res = Enumerable.Range(0, 1000000).Select(i => Sampling.SampleDiscDist(_rnd, DistType.Geometric, 62)).Average();
        Assert.AreEqual(res, 62, 0.1);
    }

    [Test]
    public void TestGammaRecoversMeanAndShape([Values(1.5, 4.0, 13.6)] double shape)
    {
        var samples = Enumerable.Range(0, 1000000)
            .Select(i => (double) Sampling.SampleDiscDist(_rnd, DistType.Gamma, 62, shape)).ToList();
        // Shape sets the dispersion only: mean must come out at RateMean regardless of it, and
        // the CV at 1/sqrt(shape) -- the property that lets one shape serve many cancer types.
        Assert.AreEqual(62, samples.Mean(), 0.1);
        Assert.AreEqual(1 / Math.Sqrt(shape), samples.StandardDeviation() / samples.Mean(), 0.01);
    }

    [Test]
    public void TestGammaAtShapeOneMatchesGeometric()
    {
        // shape 1 is the exponential, so it should reproduce the geometric SimChA sampled before.
        var gamma = Enumerable.Range(0, 200000)
            .Select(i => (double) Sampling.SampleDiscDist(_rnd, DistType.Gamma, 62, 1)).ToList();
        Assert.AreEqual(62, gamma.Mean(), 0.5);
        Assert.AreEqual(1.0, gamma.StandardDeviation() / gamma.Mean(), 0.02);
    }

    [Test]
    public void TestGammaCountsAreUsable()
    {
        // A tiny mean must still yield at least one event rather than a zero-event sample.
        var counts = Enumerable.Range(0, 10000)
            .Select(i => Sampling.SampleDiscDist(_rnd, DistType.Gamma, 0.5, 1.5)).ToList();
        Assert.GreaterOrEqual(counts.Min(), 1);
    }

    [Test]
    public void TestGammaIsSeedDeterministic()
    {
        var first = Enumerable.Range(0, 100)
            .Select(i => Sampling.SampleDiscDist(new Random(7), DistType.Gamma, 62, 1.5)).ToList();
        var second = Enumerable.Range(0, 100)
            .Select(i => Sampling.SampleDiscDist(new Random(7), DistType.Gamma, 62, 1.5)).ToList();
        Assert.AreEqual(first, second);
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
        Assert.AreEqual(119319, Sampling.GetNormSeg(_rnd, 1_000_000, 0.1));
    }

    [Test]
    public void TestDirichlet()
    {
        var mixture = Sampling.CreateRandomMixture(_rnd, [1000.0, 1.0, 0.001]);
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
            int index = rnd.PickRndIndex(probs);
            Assert.GreaterOrEqual(index, 0);
            Assert.Less(index, probs.Count);
            Assert.AreNotEqual(index, 2);
            Assert.AreNotEqual(index, 3);
        }
    }
}
