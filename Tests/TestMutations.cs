using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SimChA.DataTypes;
using SimChA.Simulation;

namespace Tests;

public class TestMutations
{
    [SetUp]
    public void Setup()
    {
    }

    [Test]
    public void TestDeletion()
    {
        var cRegion = ReferenceGenome.GetGenotype(true)[0];
        var regions = new List<Region> { cRegion };
        // Start and end within a region
        regions = RegionOps.DeleteRange(regions, 1000, 2000);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(cRegion.Length - 1000, Chromosome.Length(regions));
        // Start and end within a region (out of two)
        regions = RegionOps.DeleteRange(regions, 2000, 4000);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(cRegion.Length - 3000, Chromosome.Length(regions));
        // Two neighbouring regions
        regions = RegionOps.DeleteRange(regions, 1500, 2500);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(cRegion.Length - 4000, Chromosome.Length(regions));
        // Cut region out
        regions = RegionOps.DeleteRange(regions, 500, 2500);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(cRegion.Length - 6000, Chromosome.Length(regions));
        // Remove region from front
        regions = RegionOps.DeleteRange(regions, 0, 500);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(cRegion.Length - 6500, Chromosome.Length(regions));
        // oversized range selection
        regions = RegionOps.DeleteRange(regions, -1000, cRegion.Length);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(0, Chromosome.Length(regions));
    }

    [Test]
    public void TestDeletionUnitLength()
    {
        var cRegion = ReferenceGenome.GetGenotype(true)[0];
        var regions = new List<Region>
        {
            cRegion with { Start = 0, End = 1 }, 
            cRegion with { Start = 1, End = 2 }, 
            cRegion with { Start = 2, End = 3 }, 
            cRegion with { Start = 3, End = 4 }, 
        };
        var res = new List<Region>
        {
            cRegion with { Start = 0, End = 1 }, 
            cRegion with { Start = 3, End = 4 }
        };
        Assert.AreEqual(res, RegionOps.DeleteRange(regions, 1, 3));
    }

    [Test]
    public void TestCopy()
    {
        var cRegion = ReferenceGenome.GetGenotype(true)[0];
        var regions = new List<Region> { cRegion };
        regions = RegionOps.DeleteRange(regions, 1000, 2000);
        regions = RegionOps.DeleteRange(regions, 4000, 5000);
        var regCopy = RegionOps.CopyRange(regions, 500, 3500);
        Console.WriteLine(Chromosome.ToString(regCopy));
    }
    
    [Test]
    public void TestCopyUnitLength()
    {
        var cRegion = ReferenceGenome.GetGenotype(true)[0];
        var regions = new List<Region>
        {
            cRegion with { Start = 0, End = 1 }, 
            cRegion with { Start = 1, End = 2 }, 
            cRegion with { Start = 2, End = 3 }, 
            cRegion with { Start = 3, End = 4 }, 
        };
        var res = new List<Region>
        {
            cRegion with { Start = 1, End = 2 }, 
            cRegion with { Start = 2, End = 3 }
        };
        Assert.AreEqual(res, RegionOps.CopyRange(regions, 1, 3));
    }
}