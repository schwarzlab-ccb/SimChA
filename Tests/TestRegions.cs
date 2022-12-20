using System;
using System.Collections.Generic;
using NUnit.Framework;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.Simulation;

namespace Tests;

[TestFixture]
public class TestRegions
{
    private Region _cRegion;

    [SetUp]
    public void Setup()
    {
        _cRegion = ReferenceGenome.GetGenotype(true)[0];
    }

    [Test]
    public void TestDeletion()
    {
        var regions = new List<Region> { _cRegion };
        // Start and end within a region
        regions = RegionOps.DeleteRange(regions, 1000, 2000);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(_cRegion.Length - 1000, Chromosome.Length(regions));
        // Start and end within a region (out of two)
        regions = RegionOps.DeleteRange(regions, 2000, 4000);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(_cRegion.Length - 3000, Chromosome.Length(regions));
        // Two neighbouring regions
        regions = RegionOps.DeleteRange(regions, 1500, 2500);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(_cRegion.Length - 4000, Chromosome.Length(regions));
        // Cut region out
        regions = RegionOps.DeleteRange(regions, 500, 2500);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(_cRegion.Length - 6000, Chromosome.Length(regions));
        // Remove region from front
        regions = RegionOps.DeleteRange(regions, 0, 500);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(_cRegion.Length - 6500, Chromosome.Length(regions));
        // oversized range selection
        regions = RegionOps.DeleteRange(regions, -1000, _cRegion.Length);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(0, Chromosome.Length(regions));
    }

    [Test]
    public void TestDeletionUnitLength()
    {
        var regions = new List<Region>
        {
            _cRegion with { Start = 0, End = 1 },
            _cRegion with { Start = 1, End = 2 },
            _cRegion with { Start = 2, End = 3 },
            _cRegion with { Start = 3, End = 4 },
        };
        var res = new List<Region>
        {
            _cRegion with { Start = 0, End = 1 },
            _cRegion with { Start = 3, End = 4 }
        };
        Assert.AreEqual(res, RegionOps.DeleteRange(regions, 1, 3));
    }

    [Test]
    public void TestCopy()
    {
        var regions = new List<Region> { _cRegion };
        // Start and end within a region
        var regCopy = RegionOps.CopyRange(regions, 500, 3500);
        Console.WriteLine(Chromosome.ToString(regCopy));
        Assert.AreEqual(3000, Chromosome.Length(regCopy));
        // Copy across regions
        regions = RegionOps.DeleteRange(regions, 1000, 2000);
        regions = RegionOps.DeleteRange(regions, 4000, 5000);
        regCopy = RegionOps.CopyRange(regions, 500, 3500);
        Console.WriteLine(Chromosome.ToString(regCopy));
        Assert.AreEqual(3000, Chromosome.Length(regCopy));
    }

    [Test]
    public void TestCopyUnitLength()
    {
        var regions = new List<Region>
        {
            _cRegion with { Start = 0, End = 1 },
            _cRegion with { Start = 1, End = 2 },
            _cRegion with { Start = 2, End = 3 },
            _cRegion with { Start = 3, End = 4 },
        };
        var res = new List<Region>
        {
            _cRegion with { Start = 1, End = 2 },
            _cRegion with { Start = 2, End = 3 }
        };
        Assert.AreEqual(res, RegionOps.CopyRange(regions, 1, 3));
    }
    
    [Test]
    public void TestSplit()
    {
        var regions = new List<Region> { _cRegion };
        var (before, after) = RegionOps.SplitRegions(regions, 2000);
        Console.WriteLine(Chromosome.ToString(before));
        Assert.AreEqual(2000, Chromosome.Length(before));
        Console.WriteLine(Chromosome.ToString(after));
        Assert.AreEqual(_cRegion.Length - 2000, Chromosome.Length(after));
    }
    
    [Test]
    public void TestSplitUnitLength()
    {
        var regions = new List<Region>
        {
            _cRegion with { Start = 0, End = 2 },
        };
        var (before, after) = RegionOps.SplitRegions(regions, 1);
        Assert.AreEqual(1, Chromosome.Length(after));
    }
    
    [Test]
    public void TestInverse()
    {
        var regions = new List<Region> { _cRegion };
        regions = RegionOps.DeleteRange(regions, 1000, 2000);
        regions = RegionOps.InvertRegions(regions);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(_cRegion.Length - 1000, Chromosome.Length(regions));
    }
    
    [Test]
    public void TestConcat()
    {
        var regions = new List<Region> { _cRegion };
        regions = RegionOps.ConcatRegions(regions, regions);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(_cRegion.Length * 2, Chromosome.Length(regions));
        regions = RegionOps.ConcatRegions(new[] {regions, regions });
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(_cRegion.Length * 4, Chromosome.Length(regions));
    }
}