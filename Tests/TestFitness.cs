// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using System.Collections.Generic;
using NUnit.Framework;
using SimChA.Computation;
using SimChA.DataTypes;

namespace Tests;

[TestFixture]
public class TestFitness
{
    private static Gene MakeGene(ChrNo chrNo, double deltaFitness)
        => new("G" + chrNo, new Region(0, 50, new ChrID(chrNo, false)), deltaFitness);

    [Test]
    public void TestEssTerm()
    {
        Assert.AreEqual(0, Fitness.EssTerm(new List<(Gene, int)>()));
        
        var testNoEffect = new List<(Gene, int)> { (MakeGene(ChrNo.chr1, 0), 0) };
        Assert.AreEqual(0, Fitness.EssTerm(testNoEffect));
        
        var testMissing = new List<(Gene, int)> { (MakeGene(ChrNo.chr1, 0.1), 0) };
        Assert.AreEqual(0.1, Fitness.EssTerm(testMissing));

        var testHaplosufficient = new List<(Gene, int)> { (MakeGene(ChrNo.chr1, 0.1), 1) };
        Assert.AreEqual(0, Fitness.EssTerm(testHaplosufficient));
        
        var testList = new List<(Gene, int)> { (MakeGene(ChrNo.chr1, 0.1), 0), (MakeGene(ChrNo.chr2, 0.2), 0) };
        Assert.AreEqual(0.1 + 0.2, Fitness.EssTerm(testList));
    }

    [Test]
    public void TestTsgOgTerm()
    {
        // TODO: @Felix
    }
    
    [Test]
    public void TestStressTerm()
    {
        // TODO: @Felix
    }
    
    [Test]
    public void TestCNCalulation()
    {
        // TODO: @Felix
    }
}