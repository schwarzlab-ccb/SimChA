using System;
using System.Collections.Generic;
using NUnit.Framework;
using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Simulation;

namespace Tests;

[TestFixture]
public class TestKaryotype
{
    private Karyotype _kar;
    private Random _rnd;
    
    [SetUp]
    public void Setup()
    {
        _kar = new Karyotype(false);
        _rnd = new Random(0);
    }

    [Test]
    public void TestDeletion()
    {
        _kar.ApplyAberration(_rnd, AberrationEnum.ChromDeletion, new BaseAbbP(1f));
        Assert.AreEqual(45, _kar.ContigCount);
    }
    
    [Test]
    public void TestDuplication()
    {
        _kar.ApplyAberration(_rnd, AberrationEnum.ChromDuplication, new BaseAbbP(1f));
        Assert.AreEqual(47, _kar.ContigCount);
    }
    
    [Test]
    public void TestBFB()
    {
        _kar.ApplyAberration(_rnd, AberrationEnum.BreakageFusionBridge, new FractionAbbP(1f, .1f));
        Assert.AreEqual(46, _kar.ContigCount);
    }
    
    [Test]
    public void TestClean()
    {
        for (int i = 0; i < 46; i++)
        {
            _kar.ApplyAberration(_rnd, AberrationEnum.ChromDeletion, new BaseAbbP(1f));
        }
        Assert.AreEqual("[]", _kar.ToString());
    }

    [Test]
    public void TestGetPresentGenes()
    {
        _kar = new Karyotype(false);
        var tsgOgLists = new Dictionary<ChrNo, List<Gene>>();
        foreach(var chrNo in (ChrNo[]) Enum.GetValues(typeof(ChrNo)))
        {
            tsgOgLists.Add(chrNo, new List<Gene>{new Gene("T" + chrNo.ToString(), 
                new Region(0, 50, new ChrID(chrNo, false)), (float) _rnd.NextDouble())});
        }
        var tsgOgsPresent = _kar.GetPresentGenes(tsgOgLists);
        Assert.AreEqual(_kar.ContigCount, tsgOgsPresent.Count);
        _kar.ApplyAberration(_rnd, AberrationEnum.ChromDeletion, new BaseAbbP(1f));
        _kar.ApplyAberration(_rnd, AberrationEnum.ChromDeletion, new BaseAbbP(1f));
        tsgOgsPresent = _kar.GetPresentGenes(tsgOgLists);
        Assert.AreEqual(_kar.ContigCount, tsgOgsPresent.Count);
    }
}