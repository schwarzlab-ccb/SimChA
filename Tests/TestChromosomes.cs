// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using System.Collections.Generic;
using NUnit.Framework;
using SimChA.DataTypes;

namespace Tests;

[TestFixture]
public class TestChromosomes
{
    private Chromosome _chr1;
    private Chromosome _chrX;
    
    [SetUp]
    public void Setup()
    {
        _chr1 = new Chromosome(ReferenceGenome.GetRegion(ChromNum.chr1));
        _chrX = new Chromosome(ReferenceGenome.GetRegion(ChromNum.chrX, false));
    }

    [Test]
    public void TestSplit()
    {
        int remainderLen = _chr1.Length() - 1000;
        var rest = _chr1.Split(1000, true);
        Assert.AreEqual(1000, _chr1.Length());
        Assert.AreEqual(remainderLen, rest.Length());
    }

    [Test]
    public void TestJoin()
    {
        int combinedLen = _chr1.Length() + _chrX.Length(); 
        _chr1.Join(_chrX, true);
        Assert.AreEqual(combinedLen, _chr1.Length());
    }

    [Test]
    public void TestInversion()
    {
        int length = _chr1.Length(); 
        _chr1.InvertRange(length / 4, length * 3 / 4);
        Assert.AreEqual(length, _chr1.Length());
    }

    [Test]
    public void TestReplication()
    {
        int length = _chr1.Length() + 900;
        _chr1.DuplicateRange(100, 1000);
        Assert.AreEqual(length, _chr1.Length());
    }
    
    [Test]
    public void TestBridge()
    {
        int length = (_chr1.Length() - 1000) * 2;
        _chr1.Bridge(1000, true);
        Assert.AreEqual(length, _chr1.Length());
    }

    [Test]
    public void TestScatterAndGather()
    {
        int length = _chr1.Length();
        _chr1.ScatterAndGather(new List<int>{1000, 2000, 3000}, 4);
        Assert.AreEqual(length, _chr1.Length());
    }
}