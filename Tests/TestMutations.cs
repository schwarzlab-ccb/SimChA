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
        var chr1 = new Chromosome(ReferenceGenome.Genotype[0]);
        // Start and end within a region
        chr1.DeleteRange(1000, 2000);
        Console.WriteLine(chr1);
        Assert.AreEqual(ReferenceGenome.Genotype[0].Length - 1000, chr1.Length);
        // Start and end within a region (out of two)
        chr1.DeleteRange(2000, 4000);
        Console.WriteLine(chr1);
        Assert.AreEqual(ReferenceGenome.Genotype[0].Length - 3000, chr1.Length);
        // Two neighbouring regions
        chr1.DeleteRange(1500, 2500);
        Console.WriteLine(chr1);
        Assert.AreEqual(ReferenceGenome.Genotype[0].Length - 4000, chr1.Length);
        // Cut region out
        chr1.DeleteRange(500, 2500);
        Console.WriteLine(chr1);
        Assert.AreEqual(ReferenceGenome.Genotype[0].Length - 6000, chr1.Length);
        // Remove region from front
        chr1.DeleteRange(0, 500);
        Console.WriteLine(chr1);
        Assert.AreEqual(ReferenceGenome.Genotype[0].Length - 6500, chr1.Length);
        // oversized range selection
        chr1.DeleteRange(-1000, ReferenceGenome.Genotype[0].Length);
        Console.WriteLine(chr1);
        Assert.AreEqual(0, chr1.Length);
    }
}