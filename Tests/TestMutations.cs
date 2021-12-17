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
        var regions = new List<Region> {cRegion};
        // Start and end within a region
        regions = ChrMutations.DeleteRange(regions, 1000, 2000);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(cRegion.Length - 1000, Chromosome.Length(regions));
        // Start and end within a region (out of two)
        regions = ChrMutations.DeleteRange(regions, 2000, 4000);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(cRegion.Length - 3000, Chromosome.Length(regions));
        // Two neighbouring regions
        regions = ChrMutations.DeleteRange(regions, 1500, 2500);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(cRegion.Length - 4000, Chromosome.Length(regions));
        // Cut region out
        regions = ChrMutations.DeleteRange(regions, 500, 2500);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(cRegion.Length - 6000, Chromosome.Length(regions));
        // Remove region from front
        regions = ChrMutations.DeleteRange(regions, 0, 500);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(cRegion.Length - 6500, Chromosome.Length(regions));
        // oversized range selection
        regions = ChrMutations.DeleteRange(regions, -1000, cRegion.Length);
        Console.WriteLine(Chromosome.ToString(regions));
        Assert.AreEqual(0, Chromosome.Length(regions));
    }
    
    [Test]
    public void TestCopy()
    {
        var cRegion = ReferenceGenome.GetGenotype(true)[0];
        var regions = new List<Region> {cRegion};
        regions = ChrMutations.DeleteRange(regions, 1000, 2000);
        regions = ChrMutations.DeleteRange(regions, 4000, 5000);
        var regCopy = ChrMutations.CopyRange(regions, 500, 3500);
        Console.WriteLine(Chromosome.ToString(regCopy));
        
    }
}