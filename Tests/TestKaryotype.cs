using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;
using SimChA.Simulation;

namespace Tests;

[TestFixture]
public class TestKaryotype
{
    private Karyotype _kar;
    private Random _rnd;
    private CNEventPars _del;
    private CNEventPars _dup;
    
    private int TEST_FRAC = 1000;
    
    public static void ApplyRandomEvent(Random rnd, Karyotype kar, CNEventPars cnEventPars)
    {
        var eventData = Sampling.GenerateCNEventData(rnd, kar, cnEventPars);
        eventData.ApplyEvent(kar);
    }
    
    [SetUp]
    public void Setup()
    {
        _kar = new Karyotype(false);
        _rnd = new Random(0);
        _del = new CNEventPars(CNEventType.ChromDeletion);
        _dup = new CNEventPars(CNEventType.ChromDuplication);
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
        var gluedRegions = RegionOps.GlueNeighbours(_kar.FindRegionsOfChr(ChrNo.chr1).ToList());
        Assert.AreEqual(3, gluedRegions.Count(r => r.ChrID.Parent));
        Assert.AreEqual(1, gluedRegions.Count(r => !r.Forward));
        Assert.AreEqual(TEST_FRAC, gluedRegions.First(r => !r.Forward).Start);
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
        long contig0Len = _kar.ContigLen(0);
        long chrLen = RegionOps.GetLength(_kar.FindRegionsOfChr(ChrNo.chr1).ToList());
        
        _kar.ApplyTranslocation(0, 1, TEST_FRAC, TEST_FRAC, false);
        Assert.AreEqual(contig0Len, _kar.ContigLen(1));
        Assert.AreEqual(chrLen, RegionOps.GetLength(_kar.FindRegionsOfChr(ChrNo.chr1).ToList()));

        _kar.ApplyTranslocation(0, 1, TEST_FRAC, TEST_FRAC, true);
        Assert.AreEqual(chrLen, RegionOps.GetLength(_kar.FindRegionsOfChr(ChrNo.chr1).ToList()));
        Assert.AreEqual(contig0Len, _kar.ContigLen(0));
        Console.WriteLine(_kar);
    }
    
    [Test]
    public void TestApplyCNEvent()
    {
        var pars = new Dictionary<string, double> { ["Size"] = 1_000_000 };
        Assert.DoesNotThrow(() => { ApplyRandomEvent(_rnd, _kar, new CNEventPars(CNEventType.ChromDeletion, 1.0)); });
        Assert.DoesNotThrow(() => { ApplyRandomEvent(_rnd, _kar, new CNEventPars(CNEventType.ChromDuplication, 1.0)); });
        Assert.DoesNotThrow(() => { ApplyRandomEvent(_rnd, _kar, new CNEventPars(CNEventType.InternalDeletion, 1.0, pars)); });
        Assert.DoesNotThrow(() => { ApplyRandomEvent(_rnd, _kar, new CNEventPars(CNEventType.InternalDuplication, 1.0, pars)); });
        Assert.DoesNotThrow(() => { ApplyRandomEvent(_rnd, _kar, new CNEventPars(CNEventType.InternalInversion, 1.0, pars)); });
        Assert.DoesNotThrow(() => { ApplyRandomEvent(_rnd, _kar, new CNEventPars(CNEventType.InvertedDuplication, 1.0, pars)); });
        Assert.DoesNotThrow(() => { ApplyRandomEvent(_rnd, _kar, new CNEventPars(CNEventType.BreakageFusionBridge, 1.0, pars)); });
        Assert.DoesNotThrow(() => { ApplyRandomEvent(_rnd, _kar, new CNEventPars(CNEventType.TailDeletion, 1.0, pars)); });
        Assert.DoesNotThrow(() => { ApplyRandomEvent(_rnd, _kar, new CNEventPars(CNEventType.Translocation, 1.0)); });
        Assert.DoesNotThrow(() => { ApplyRandomEvent(_rnd, _kar, new CNEventPars(CNEventType.Chromoplexy, 1.0)); });
        Assert.DoesNotThrow(() => { ApplyRandomEvent(_rnd, _kar, new CNEventPars(CNEventType.Chromothripsis, 1.0)); });
        Assert.DoesNotThrow(() => { ApplyRandomEvent(_rnd, _kar, new CNEventPars(CNEventType.Pyrgo, 1.0)); });
        Assert.DoesNotThrow(() => { ApplyRandomEvent(_rnd, _kar, new CNEventPars(CNEventType.Rigma, 1.0)); });
        Assert.DoesNotThrow(() => { ApplyRandomEvent(_rnd, _kar, new CNEventPars(CNEventType.TIChain, 1.0)); });
        Assert.DoesNotThrow(() => { ApplyRandomEvent(_rnd, _kar, new CNEventPars(CNEventType.TIBridge, 1.0)); });
        Assert.DoesNotThrow(() => { ApplyRandomEvent(_rnd, _kar, new CNEventPars(CNEventType.TICycle, 1.0)); });
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
            new()
        };
        var sequence = new List<int> { 0, 7, 6, 5, 4, 3, 2, 1 };
        var breakpoints = new List<long> { _kar.ContigLen(2), _kar.ContigLen(1) };
        _kar.ApplyChromoplexy(ids, stops, sequence, breakpoints);
        Assert.AreEqual(46, _kar.CountContigs());
        Assert.AreEqual(HGRef.GetGenomeLen(_kar.SexXX), _kar.GenomeLen());
    }
    
    [Test]
    public void TestClean()
    {
        for (int i = 0; i < HGRef.CHR_COUNT; i++)
        {
            ApplyRandomEvent(_rnd, _kar, _del);
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
        ApplyRandomEvent(_rnd, _kar, _del);
        tsgOgsPresent = _kar.GetPresentGenes(tsgOgLists);
        Assert.AreEqual(_kar.CountContigs(), tsgOgsPresent.Count);
    }
    
    [Test]
    public void TestPyrgo()
    {
        long contigLen = _kar.ContigLen(0);
        var frags = new List<(long, long)> { (TEST_FRAC, TEST_FRAC), (TEST_FRAC, TEST_FRAC * 2) };
        _kar.ApplyPyrgo(0, frags);
        Assert.AreEqual(contigLen, _kar.ContigLen(0) - TEST_FRAC * 3);
    }
    
    [Test]
    public void TestRigma()
    {
        long contigLen = _kar.ContigLen(0);
        _kar.ApplyRigma(0, TEST_FRAC, new List<long> {TEST_FRAC, TEST_FRAC, TEST_FRAC});
        Assert.AreEqual(contigLen, _kar.ContigLen(0) + TEST_FRAC * 2);
    }

    [Test]
    public void TestTIBridge()
    {
        long contigLen = _kar.ContigLen(0);
        var frags = new List<(int, long, long, bool)>
        {
            (0, TEST_FRAC, 0, true),
            (1, TEST_FRAC, TEST_FRAC, false),
            (2, TEST_FRAC * 2, TEST_FRAC * 2, true)
        };
        _kar.ApplyTIBridge(frags);
        Assert.AreEqual(46, _kar.CountContigs());
        Assert.AreEqual(contigLen + TEST_FRAC * 3, _kar.ContigLen(0));
    }

    [Test]
    public void TestTIChain()
    {
        var frags = new List<(int, long, long, bool)>
        {
            (0, TEST_FRAC, TEST_FRAC, true),
            (1, TEST_FRAC, TEST_FRAC * 2, false),
            (2, TEST_FRAC * 2, TEST_FRAC, true)
        };
        _kar.ApplyTIChain(frags);
        Assert.AreEqual(47, _kar.CountContigs());
        Assert.AreEqual(TEST_FRAC * 4, _kar.ContigLen(46));
    }

    [Test]
    public void TestTICycle()
    {
        long contigLen = _kar.ContigLen(0);
        var frags = new List<(int, long, long, bool)>
        {
            (0, TEST_FRAC, TEST_FRAC, true),
            (1, TEST_FRAC, TEST_FRAC, false),
            (2, TEST_FRAC * 2, TEST_FRAC * 2, true)
        };
        _kar.ApplyTICycle(frags);
        Assert.AreEqual(46, _kar.CountContigs());
        Assert.AreEqual(contigLen + TEST_FRAC * 4, _kar.ContigLen(0));
    }
}
