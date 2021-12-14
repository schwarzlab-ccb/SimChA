using System;
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
        ChrMutations.DeleteRegion(chr1, 1000, 2000);
        Console.WriteLine(chr1);
        Assert.AreEqual(ReferenceGenome.Genotype[0].Length, chr1.Length + 1000);
    }
    
    [Test]
    public void TestDeletion2()
    {
        var chr1 = new Chromosome(ReferenceGenome.Genotype[0]);
        ChrMutations.DeleteRegion(chr1, 0, 1000);
        Console.WriteLine(chr1);
        Assert.AreEqual(ReferenceGenome.Genotype[0].Length, chr1.Length + 1000);
    }
}