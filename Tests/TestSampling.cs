using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SimChA.Simulation;

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
        var probs = new List<double> {0.1, 0.2, 0.3, 0.4};
        int index = Sampling.PickRandomIndex(_rnd, probs);
        Assert.GreaterOrEqual(index, 0);
        Assert.Less(index, probs.Count);
    }
}
