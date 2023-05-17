// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using System.Collections.Generic;
using NUnit.Framework;
using SimChA.DataTypes;
using SimChA.Simulation;

namespace Tests;

[TestFixture]
public class TestContig
{
    private Contig _contig1;
    private Contig _contigX;
    
    [SetUp]
    public void Setup()
    {
        HGRef.Assembly = GenomeAssembly.hg19;
        _contig1 = new Contig(HGRef.GetRegion(ChrNo.chr1));
        _contigX = new Contig(HGRef.GetRegion(ChrNo.chrX, false));
    }
    
    [Test]
    public void TestSplit()
    {
        long remainderLen = _contig1.Length() - 1000;
        var rest = _contig1.Split(1000, true);
        Assert.AreEqual(1000, _contig1.Length());
        Assert.AreEqual(remainderLen, rest.Length());
    }

    [Test]
    public void TestJoin()
    {
        long combinedLen = _contig1.Length() + _contigX.Length(); 
        _contig1.Join(_contigX);
        Assert.AreEqual(combinedLen, _contig1.Length());
    }

    [Test]
    public void TestInversion()
    {
        long length = _contig1.Length(); 
        _contig1.InvertRange(length / 4, length * 3 / 4);
        Assert.AreEqual(length, _contig1.Length());
    }

    [Test]
    public void TestReplication()
    {
        long length = _contig1.Length() + 900;
        _contig1.DuplicateRange(100, 1000);
        Assert.AreEqual(length, _contig1.Length());
    }
    
    [Test]
    public void TestBridgeFront()
    {
        long length = (_contig1.Length() - 1000) * 2;
        _contig1.Bridge(1000, true);
        Assert.AreEqual(length, _contig1.Length());
        _contig1.Bridge(1000, false);
        Assert.AreEqual(1000*2, _contig1.Length());
    }

    [Test]
    public void TestScatterAndGather()
    {
        long length = _contig1.Length();
        _contig1.ScatterAndGather(new List<long>{1000, 2000, 3000}, new List<int>{3, 1, 2, 0});
        Assert.AreEqual(length, _contig1.Length());
    }

    [Test]
    public void TestGetRandomRegion()
    {
        var res = _contig1.GetSubContig(1, _contig1.Length() - 1);
        Assert.LessOrEqual(_contig1.Length() - 2, res.Length());
    }

    [Test]
    public void TestInsertContig()
    {
        var copyOfContig1 = new Contig(_contig1);
        copyOfContig1.InsertContig(_contigX, _contig1.Length() / 2);
        Assert.AreEqual(_contig1.Length() + _contigX.Length(), copyOfContig1.Length());
    }
}