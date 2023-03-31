using System;
using System.Collections.Generic;
using NUnit.Framework;
using SimChA.Simulation;

namespace Tests;

[TestFixture]

public class TestSampling
{
    [Test]
    public void TestGetExpSeg()
    {
        var rnd = new Random(0);
        Assert.AreEqual(31987, Sampling.GetExpSeg(rnd, 1_000_000, 0.1));
        Assert.AreEqual(202, Sampling.GetExpSeg(rnd, 1_000_000, 1_000));
    }
    
    [Test]
    public void TestGetNormSeg()
    {
        var rnd = new Random(0);
        Assert.AreEqual(127097, Sampling.GetNormSeg(rnd, 1_000_000, 0.1));
    }
}
