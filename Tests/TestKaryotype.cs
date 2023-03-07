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
    private CNEventP _del;
    private CNEventP _dup;
    
    private int TEST_FRAC = 1000;
    
    [SetUp]
    public void Setup()
    {
        _kar = new Karyotype(false);
        _rnd = new Random(0);
        _del = new CNEventP(CNEventType.ChromDeletion, 1f);
        _dup = new CNEventP(CNEventType.ChromDuplication, 1f);
    }

    // Test for each AberrationEnum value
    [Test]
    public void TestWGD()
    {
        _kar.ApplyWGD();
        Assert.AreEqual(92, _kar.CountContigs());
    }

    [Test]
    public void TestContigDeletion()
    {
        _kar.ApplyContigDeletion(0);
        Assert.AreEqual(45, _kar.CountContigs());
        
        _kar.ApplyContigDeletion(45);
        Assert.AreEqual(44, _kar.CountContigs());
    }
    
    [Test]
    public void TestContigDuplication()
    {
        _kar.ApplyContigDuplication(0);
        Assert.AreEqual(47, _kar.CountContigs());
        
        _kar.ApplyContigDuplication(46);
        Assert.AreEqual(48, _kar.CountContigs());
    }

    [Test]
    public void TestInternalDeletion()
    {
        long len = _kar.ContigLen(0);
        _kar.ApplyInternalDeletion(0, TEST_FRAC, 2 * TEST_FRAC);
        Assert.AreEqual(len - TEST_FRAC, _kar.ContigLen(0));
        _kar.ApplyInternalDeletion(0, TEST_FRAC, 2 * TEST_FRAC);
        Assert.AreEqual(len - 2 * TEST_FRAC, _kar.ContigLen(0));
    }
    
    [Test]
    public void TestInternalDuplication()
    {
        long len = _kar.ContigLen(0);
        _kar.ApplyInternalDuplication(0, TEST_FRAC, 2 * TEST_FRAC);
        Assert.AreEqual(len + TEST_FRAC, _kar.ContigLen(0));
        _kar.ApplyInternalDuplication(0, TEST_FRAC, 2 * TEST_FRAC);
        Assert.AreEqual(len + 2 * TEST_FRAC, _kar.ContigLen(0));
    }

    [Test]
    public void TestInternalInversion()
    {
        long len = _kar.ContigLen(0);
        _kar.ApplyInternalInversion(0, TEST_FRAC, 2 * TEST_FRAC);
        Assert.AreEqual(len, _kar.ContigLen(0));
        var regions = _kar.FindRegionsOfChr(ChrNo.chr1).ToList();
        Assert.AreEqual(3, regions.Count(r => r.ChrID.Parent));
        Assert.AreEqual(1, regions.Count(r => !r.Forward));
        Assert.AreEqual(TEST_FRAC, regions.First(r => !r.Forward).Start);
    }
    
    [Test]
    public void TestInvertedDuplication()
    {
        long len = _kar.ContigLen(0);
        _kar.ApplyInvertedDuplication(0, TEST_FRAC, 2 * TEST_FRAC);
        Assert.AreEqual(len + TEST_FRAC, _kar.ContigLen(0));
        var regions = _kar.FindRegionsOfChr(ChrNo.chr1).ToList();
        Assert.AreEqual(3, regions.Count(r => r.ChrID.Parent));
        Assert.AreEqual(1, regions.Count(r => !r.Forward));
        Assert.AreEqual(TEST_FRAC, regions.First(r => !r.Forward).Start);
    }
    
    [Test]
    public void TestTailDeletion()
    {
        long len = _kar.ContigLen(0);
        _kar.ApplyTailDeletion(0, TEST_FRAC, true);
        Assert.AreEqual(len - TEST_FRAC, _kar.ContigLen(0));
        _kar.ApplyTailDeletion(0, TEST_FRAC, false);
        Assert.AreEqual(len - 2 * TEST_FRAC, _kar.ContigLen(0));
    }

    [Test]
    public void TestBFB()
    {
        long len = _kar.ContigLen(0);
        _kar.ApplyBFB(0, TEST_FRAC, true);
        long newLen = _kar.ContigLen(0);
        Assert.AreEqual((len - TEST_FRAC) * 2, newLen);
    }
    
    [Test]
    public void TestTranslocation()
    {
        long contigLen = _kar.ContigLen(0);
        long chrLen = RegionOps.GetLength(_kar.FindRegionsOfChr(ChrNo.chr1).ToList());
        _kar.ApplyTranslocation(0, 1, TEST_FRAC, TEST_FRAC);
        Assert.AreEqual(contigLen, _kar.ContigLen(1));
        Assert.AreEqual(chrLen, RegionOps.GetLength(_kar.FindRegionsOfChr(ChrNo.chr1).ToList()));
    }

    [Test]
    public void TestChromothripsis()
    {
        long contigLen = _kar.ContigLen(0);
        var stops = new List<long> { TEST_FRAC, TEST_FRAC * 2, TEST_FRAC * 3 };
        var selection = new List<int> { 3, 1}; // Keep only a TEST_FRAC chunk and the tail
        _kar.ApplyChromothripsis(0, stops, selection);
        Assert.AreEqual(contigLen - TEST_FRAC * 2, _kar.ContigLen(0));
    }
    
    [Test]
    public void TestClean()
    {
        for (int i = 0; i < HGRef.CHR_COUNT; i++)
        {
            _kar.ApplyAberration(_rnd, _del);
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
        Assert.AreEqual(_kar.CountContigs(), tsgOgsPresent.Count);
        
        // Removes a gene and a contig at the same time
        _kar.ApplyAberration(_rnd, _del);
        tsgOgsPresent = _kar.GetPresentGenes(tsgOgLists);
        Assert.AreEqual(_kar.CountContigs(), tsgOgsPresent.Count);
    }
}