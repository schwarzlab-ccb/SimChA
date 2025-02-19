// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SimChA.Data;
using SimChA.IO;

namespace Tests;

[TestFixture]
public class TestContig
{
    private GenRef _genRef;
    private Contig _contig1;
    private Contig _contigX;
    
    [SetUp]
    public void Setup()
    {
        _genRef = FileIO.ReadGenRef("./../../../../data/hg19");
        _contig1 = new Contig(_genRef.GetGenotype(SexType.Female).First());
        _contigX = new Contig(_genRef.GetGenotype(SexType.Female).Last());
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
        for (int i = 0; i < 10; i++)
        {
            _contig1.InvertRange(i * 1000000, i * 2000000);
        }

        Assert.AreEqual(length, _contig1.Length());
    }

    [Test]
    public void TestReplication()
    {
        long length = _contig1.Length() + 10 * 900;
        for (int i = 0; i < 10; i++)
        {
            _contig1.DuplicateRange(i*1000 + 100, i*1000 + 1000);
        }
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

    [Test]
    public void TestGetCentromeres()
    {
        var cents = _contig1.GetCentromeres(_genRef.Centromeres);
        Assert.AreEqual(1, cents.Count);
        Assert.AreEqual(121500000, cents[0].start);;
        Assert.AreEqual(128900000, cents[0].end);
        
        _contig1.Bridge(1000, true);
        cents = _contig1.GetCentromeres(_genRef.Centromeres);
        Assert.AreEqual(2, cents.Count);
    }
}