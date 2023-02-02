using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SimChA.Computation;
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

    // Test for each AberrationEnum value
    [Test]
    public void TestWGD()
    {
        _kar.ApplyWGD();
        Assert.AreEqual(92, _kar.ContigCount);
    }

    [Test]
    public void TestContigDeletion()
    {
        _kar.ApplyContigDeletion(0);
        Assert.AreEqual(45, _kar.ContigCount);
        
        _kar.ApplyContigDeletion(45);
        Assert.AreEqual(44, _kar.ContigCount);
    }
    
    [Test]
    public void TestContigDuplication()
    {
        _kar.ApplyContigDuplication(0);
        Assert.AreEqual(47, _kar.ContigCount);
        
        _kar.ApplyContigDuplication(46);
        Assert.AreEqual(48, _kar.ContigCount);
    }

    [Test]
    public void TestInternalDeletion()
    {
        long len = RegionOps.GetLength(_kar.FindRegionsOfChr(ChrNo.chr1).ToList());
        _kar.ApplyInternalDeletion(0, 1000, 2000);
        Assert.AreEqual(len - 1000, RegionOps.GetLength(_kar.FindRegionsOfChr(ChrNo.chr1).ToList()));
        _kar.ApplyInternalDeletion(0, 1000, 2000);
        Assert.AreEqual(len - 2000, RegionOps.GetLength(_kar.FindRegionsOfChr(ChrNo.chr1).ToList()));
    }
    
    [Test]
    public void TestInternalDuplication()
    {
        long len = RegionOps.GetLength(_kar.FindRegionsOfChr(ChrNo.chr1).ToList());
        _kar.ApplyInternalDuplication(0, 1000, 2000);
        Assert.AreEqual(len + 1000, RegionOps.GetLength(_kar.FindRegionsOfChr(ChrNo.chr1).ToList()));
        _kar.ApplyInternalDuplication(0, 1000, 2000);
        Assert.AreEqual(len + 2000, RegionOps.GetLength(_kar.FindRegionsOfChr(ChrNo.chr1).ToList()));
    }
    
    [Test]
    public void TestInternalInversion()
    {
        long len = RegionOps.GetLength(_kar.FindRegionsOfChr(ChrNo.chr1).ToList());
        _kar.ApplyInternalInversion(0, 1000, 2000);
        Assert.AreEqual(len, RegionOps.GetLength(_kar.FindRegionsOfChr(ChrNo.chr1).ToList()));
        _kar.ApplyInternalInversion(0, 1000, 2000);
        Assert.AreEqual(len, RegionOps.GetLength(_kar.FindRegionsOfChr(ChrNo.chr1).ToList()));
    }
    
    [Test]
    public void TestInternalTranslocation()
    {
        long len = RegionOps.GetLength(_kar.FindRegionsOfChr(ChrNo.chr1).ToList());
        _kar.ApplyTailDeletion(0, 1000, true);
        Assert.AreEqual(len - 1000, RegionOps.GetLength(_kar.FindRegionsOfChr(ChrNo.chr1).ToList()));
        _kar.ApplyTailDeletion(0, 1000, false);
        Assert.AreEqual(len - 2000, RegionOps.GetLength(_kar.FindRegionsOfChr(ChrNo.chr1).ToList()));
    }

    //
    // [Test]
    // public void TestDeletion()
    // {
    //     _kar.ApplyAberration(_rnd, AberrationEnum.ChromDeletion, new BaseAbbP(1f));
    //     Assert.AreEqual(45, _kar.ContigCount);
    // }
    //
    // [Test]
    // public void TestDuplication()
    // {
    //     _kar.ApplyAberration(_rnd, AberrationEnum.ChromDuplication, new BaseAbbP(1f));
    //     Assert.AreEqual(47, _kar.ContigCount);
    // }
    //
    // [Test]
    // public void TestBFB()
    // {
    //     _kar.ApplyAberration(_rnd, AberrationEnum.BreakageFusionBridge, new FractionAbbP(1f, .1f));
    //     Assert.AreEqual(46, _kar.ContigCount);
    // }
    
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