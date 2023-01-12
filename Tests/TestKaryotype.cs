using System;
using System.Collections.Generic;
using System.Linq;
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

    private Gene MakeGene(ChrNo chrNo) 
        => new("G" + chrNo, new Region(0, 50, new ChrID(chrNo, false)), _rnd.NextDouble());

    [Test]
    public void TestGetPresentGenes()
    {
        var chrNums = Enum.GetValues(typeof(ChrNo)).Cast<ChrNo>();
        var tsgOgLists = chrNums.ToDictionary(c => c, c => new List<Gene> {MakeGene(c)});
        var tsgOgsPresent = _kar.GetPresentGenes(tsgOgLists);
        Assert.AreEqual(_kar.ContigCount, tsgOgsPresent.Count);
        
        // Removes a gene and a contig at the same time
        _kar.ApplyAberration(_rnd, AberrationEnum.ChromDeletion, new BaseAbbP(1f));
        tsgOgsPresent = _kar.GetPresentGenes(tsgOgLists);
        Assert.AreEqual(_kar.ContigCount, tsgOgsPresent.Count);
    }
}