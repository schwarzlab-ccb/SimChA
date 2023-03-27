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
        var desc = _kar.ToString();
        long contigLen = _kar.ContigLen(0);
        long chrLen = RegionOps.GetLength(_kar.FindRegionsOfChr(ChrNo.chr1).ToList());
        
        _kar.ApplyTranslocation(0, 1, TEST_FRAC, TEST_FRAC, false);
        Assert.AreEqual(contigLen, _kar.ContigLen(1));
        Assert.AreEqual(chrLen, RegionOps.GetLength(_kar.FindRegionsOfChr(ChrNo.chr1).ToList()));

        _kar.ApplyTranslocation(0, 1, TEST_FRAC, TEST_FRAC, true);
        Assert.AreEqual(contigLen, _kar.ContigLen(0));

        _kar.ApplyInternalInversion(0, 0, TEST_FRAC);
        _kar.GlueNeighbours();
        Assert.AreEqual(desc, _kar.ToString());
    }
    
    [Test]
    public void TestApplyCNEvent()
    {
        var pars = new Dictionary<string, double> { ["Mean"] = 0.1 };
        Assert.DoesNotThrow(() => { _kar.ApplyCNEvent(_rnd, new CNEventP(CNEventType.ChromDeletion, 1.0)); });
        Assert.DoesNotThrow(() => { _kar.ApplyCNEvent(_rnd, new CNEventP(CNEventType.ChromDuplication, 1.0)); });
        Assert.DoesNotThrow(() => { _kar.ApplyCNEvent(_rnd, new CNEventP(CNEventType.InternalDeletion, 1.0, pars)); });
        Assert.DoesNotThrow(() => { _kar.ApplyCNEvent(_rnd, new CNEventP(CNEventType.InternalDuplication, 1.0, pars)); });
        Assert.DoesNotThrow(() => { _kar.ApplyCNEvent(_rnd, new CNEventP(CNEventType.InternalInversion, 1.0, pars)); });
        Assert.DoesNotThrow(() => { _kar.ApplyCNEvent(_rnd, new CNEventP(CNEventType.InvertedDuplication, 1.0, pars)); });
        Assert.DoesNotThrow(() => { _kar.ApplyCNEvent(_rnd, new CNEventP(CNEventType.BreakageFusionBridge, 1.0, pars)); });
        Assert.DoesNotThrow(() => { _kar.ApplyCNEvent(_rnd, new CNEventP(CNEventType.TailDeletion, 1.0, pars)); });
        Assert.DoesNotThrow(() => { _kar.ApplyCNEvent(_rnd, new CNEventP(CNEventType.Translocation, 1.0)); });
        Assert.DoesNotThrow(() => { _kar.ApplyCNEvent(_rnd, new CNEventP(CNEventType.Chromoplexy, 1.0)); });
        Assert.DoesNotThrow(() => { _kar.ApplyCNEvent(_rnd, new CNEventP(CNEventType.Chromothripsis, 1.0)); });
        Assert.DoesNotThrow(() => { _kar.ApplyCNEvent(_rnd, new CNEventP(CNEventType.WholeGenomeDoubling, 1.0)); });
        Assert.DoesNotThrow(() => { _kar.ApplyCNEvent(_rnd, new CNEventP(CNEventType.ChainTemplatedInsertions, 1.0)); });
        Assert.DoesNotThrow(() => { _kar.ApplyCNEvent(_rnd, new CNEventP(CNEventType.BridgeTemplatedInsertions, 1.0)); });
        Assert.DoesNotThrow(() => { _kar.ApplyCNEvent(_rnd, new CNEventP(CNEventType.CycleTemplatedInsertions, 1.0)); });
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
    public void TestChromoplexy()
    {
        var ids = new List<int> { 0, 1, 2 };
        var stops = new List<List<long>>
        {
            new() { TEST_FRAC * 1, TEST_FRAC * 2 },
            new() { TEST_FRAC * 3, TEST_FRAC * 4, TEST_FRAC * 5 },
            new() { TEST_FRAC * 6, TEST_FRAC * 7 },
        };
        var sequence = new List<int> { 0, 9, 8, 7, 6, 5, 4, 3, 2, 1 };
        var breakpoints = new List<long> { _kar.ContigLen(2), _kar.ContigLen(1) };
        var result = _kar.ApplyChromoplexy(ids, stops, sequence, breakpoints);
        Assert.AreEqual("contigs:[0,1,2];fragments:10", result);
        Assert.AreEqual(46, _kar.CountContigs());
        Assert.AreEqual(HGRef.GetGenomeLen(_kar.SexXX), _kar.GenomeLen());
    }
    
    [Test]
    public void TestClean()
    {
        for (int i = 0; i < HGRef.CHR_COUNT; i++)
        {
            _kar.ApplyCNEvent(_rnd, _del);
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
        _kar.ApplyCNEvent(_rnd, _del);
        tsgOgsPresent = _kar.GetPresentGenes(tsgOgLists);
        Assert.AreEqual(_kar.CountContigs(), tsgOgsPresent.Count);
    }

    [Test]
    public void TestBridgeTemplatedInsertion()
    {
        var regions = new List<Region>();
        regions.Add(new Region(1000, 2000, new ChrID(ChrNo.chr1, false), false));
        regions.Add(new Region(3000, 5000, new ChrID(ChrNo.chr2, true), true));
        _kar.ApplyBridgeTemplatedInsertions(0, regions);
        Assert.AreEqual(_kar.ContigLen(0), HGRef.GetChromLen(ChrNo.chr1) + 3000);
    }

    [Test]
    public void TestChainTemplatedInsertion()
    {
        var regions = new List<Region>();
        regions.Add(new Region(1000, 2000, new ChrID(ChrNo.chr1, false), false));
        regions.Add(new Region(7000, 8000, new ChrID(ChrNo.chr2, true), true));
        _kar.ApplyInternalDuplication(1, 6000, 7000);
        _kar.ApplyChainTemplatedInsertions(0, regions, 1);
        Assert.AreEqual(_kar.ContigLen(0), HGRef.GetChromLen(ChrNo.chr1) + 9000);
    }

    [Test]
    public void TestCycleTemplatedInsertions()
    {
        var regions = new List<Region>();
        regions.Add(new Region(1000, 2000, new ChrID(ChrNo.chr1, false), false));
        regions.Add(new Region(7000, 8000, new ChrID(ChrNo.chr2, false), true));
        _kar.ApplyCycleTemplatedInsertions(0, regions, new Random());
        Assert.Greater(_kar.ContigLen(0), HGRef.GetChromLen(ChrNo.chr1) + 2000);
    }
}