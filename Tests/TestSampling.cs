using System;
using System.Collections.Generic;
using System.Linq;
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
    [Test]
    public void TestSampleContigsByArms()
    {
        var regions = new List<Region>
        {
            new PArm(0, 1, "chr1", true, true),
            new Centromere(1, 2, "chr1", true, true),
            new QArm(2, 3, "chr1", true, true),
            new Centromere(3, 4, "chr1", true, true),
            new QArm(4, 5, "chr1", true, true),
        };
        var contigs = new List<Contig>(){new(regions)};
        var kar = new Karyotype(contigs, new List<GenRange>(), false);
        /*for (int i = 1; i < 46; i++)
        {
            kar.ApplyContigDeletion(i);
        }*/
        _rnd = new Random(0);
        var (id, index, pArm) = Sampling.SampleContigByArms(_rnd, kar);
        Assert.AreEqual(0, id);
        // Delete the p-arm of the first chromosome
        kar.ApplyArmDeletion(0, 1, true);
        (id, index, pArm) = Sampling.SampleContigByArms(_rnd, kar);
        Assert.AreEqual(0, id);
        Assert.IsFalse(pArm);
    }
}
