using System;
using System.Collections.Generic;
using System.Linq;
using CommandLine;
using NUnit.Framework;
using SimChA.Computation;
using SimChA.EventData;
using SimChA.IO;
using SimChA.Simulation;
using SimChA.Data;

namespace Tests;

[TestFixture]
public class TestKaryotype
{
    private Karyotype _kar = null!;
    private Random _rnd = null!;
    private CNEventPars _del = null!;
    private RefGen _refGen = null!;
    private const int TEST_FRAC = 1000;
    
    public static void ApplyRandomEvent(Random rnd, Karyotype kar, CNEventPars cnEventPars)
    {
        var eventData = Sampling.GenerateCNEventData(rnd, kar, cnEventPars);
        if (eventData == null)
        {
            throw new Exception("Could not generate event data.");
        }
        eventData.ApplyEvent(kar);
    }

    // Some terminal-arm tests intentionally need a synthetic dicentric contig. Build that fixture
    // directly instead of abusing translocation, whose contract now forbids dicentric derivatives.
    private void FuseContigsForTest(int targetId, int donorId)
    {
        _kar.GetContig(targetId).Join(_kar.GetContig(donorId));
        _kar.GetContig(donorId).Clear();
    }

    
    [SetUp]
    public void Setup()
    {
        _refGen = FileIO.ReadGenRef(TestParsing.DATA_PATH, TestParsing.HG_19, TestParsing.GENE_SET);
        _kar = new Karyotype(_refGen, SexType.Male);
        _rnd = new Random(0);
        _del = new CNEventPars(CNEventType.ContigDeletion, 1);
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

    [TestCase(CNEventType.ContigDeletion)]
    [TestCase(CNEventType.ContigDuplication)]
    public void TestContigEventsCanApplyToNonChromosomeContigs(CNEventType eventType)
    {
        const int contigId = 0;
        foreach (int id in _kar.ContigIds().Where(id => id != contigId).ToList())
        {
            _kar.ApplyContigDeletion(id);
        }
        _kar.ApplyTailDeletion(contigId, 1, true);
        Assert.AreEqual(1, _kar.CountIntactTelomereEnds(contigId));
        Assert.Greater(_kar.CountCentromeres(contigId), 0);

        var eventParameters = new CNEventPars(eventType, 1);
        var eventData = Sampling.GenerateCNEventData(_rnd, _kar, eventParameters);
        Assert.NotNull(eventData);

        eventData!.ApplyEvent(_kar);
        Assert.AreEqual(
            eventType == CNEventType.ContigDeletion ? 0 : 2,
            _kar.CountContigs());
    }

    [TestCase(CNEventType.ChromDeletion)]
    [TestCase(CNEventType.ChromDuplication)]
    public void TestChromEventsRequireTwoTerminalTelomeresAndACentromere(CNEventType eventType)
    {
        const int validId = 0;
        const int oneEndedId = 1;
        foreach (int id in _kar.ContigIds().Where(id => id != validId && id != oneEndedId).ToList())
        {
            _kar.ApplyContigDeletion(id);
        }
        _kar.ApplyTailDeletion(oneEndedId, 1, true);
        Assert.AreEqual(2, _kar.CountIntactTelomereEnds(validId));
        Assert.AreEqual(1, _kar.CountIntactTelomereEnds(oneEndedId));

        var eventParameters = new CNEventPars(eventType, 1);
        for (int i = 0; i < 100; i++)
        {
            var eventData = Sampling.GenerateCNEventData(_rnd, _kar, eventParameters) as ContigEventData;
            Assert.NotNull(eventData);
            Assert.AreEqual(validId, eventData!.ContigId);
        }

        _kar.ApplyContigDeletion(oneEndedId);
        var centromere = _kar.GetCentromeres(validId).Single();
        _kar.ApplyInternalDeletion(validId, centromere.start, centromere.end);
        Assert.AreEqual(2, _kar.CountIntactTelomereEnds(validId));
        Assert.AreEqual(0, _kar.CountCentromeres(validId));
        Assert.IsNull(Sampling.GenerateCNEventData(_rnd, _kar, eventParameters));
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
        int nRegions = _kar.GetContig(0).CountRegions();
        _kar.ApplyInternalInversion(0, TEST_FRAC, 2 * TEST_FRAC);
        Assert.AreEqual(len, _kar.ContigLen(0));
        var contig = _kar.GetContig(0);
        var regions = _kar.FindChrRegions("chr1").ToList();
        Assert.AreEqual(nRegions + 2, regions.Count(r => r.Hap1));
        Assert.AreEqual(1, regions.Count(r => !r.Forward));
        Assert.AreEqual(-2 * TEST_FRAC, regions.First(r => !r.Forward).Start);
        Assert.AreEqual(TEST_FRAC, regions.First(r => !r.Forward).AbsStart);
        Assert.AreEqual(-TEST_FRAC, regions.First(r => !r.Forward).End);
        Assert.AreEqual(2 * TEST_FRAC, regions.First(r => !r.Forward).AbsEnd);
    }
    
    [Test]
    public void TestInvertedDuplication()
    {
        long len = _kar.ContigLen(0);
        _kar.ApplyInvertedDuplication(0, TEST_FRAC, 2 * TEST_FRAC);
        Assert.AreEqual(len + TEST_FRAC, _kar.ContigLen(0));
        _kar.ApplyInvertedDuplication(0, TEST_FRAC, 2 * TEST_FRAC);
        Assert.AreEqual(len + 2 * TEST_FRAC, _kar.ContigLen(0));
    }
    
    [Test]
    public void TestTailDeletion()
    {
        long len = _kar.ContigLen(0);
        _kar.ApplyTailDeletion(0, TEST_FRAC, true);
        Assert.AreEqual(len - TEST_FRAC, _kar.ContigLen(0));
        _kar.ApplyTailDeletion(0, len - 2 * TEST_FRAC, false);
        Assert.AreEqual(len - 2 * TEST_FRAC, _kar.ContigLen(0));
    }

    [Test]
    public void TestTailDuplication()
    {
        long len = _kar.ContigLen(0);
        int contigCount = _kar.CountContigs();
        _kar.ApplyTailDuplication(0, TEST_FRAC, true);
        Assert.AreEqual(contigCount, _kar.CountContigs());
        Assert.AreEqual(len + TEST_FRAC, _kar.ContigLen(0));
        _kar.ApplyTailDuplication(0, 2*TEST_FRAC, true);
        Assert.AreEqual(contigCount, _kar.CountContigs());
        Assert.AreEqual(len + 3 * TEST_FRAC, _kar.ContigLen(0));
    }

    [Test]
    public void TestTelomereEventRequiresAnIntactTerminalTelomere()
    {
        const int contigId = 0;
        Assert.AreEqual(2, _kar.CountIntactTelomereEnds(contigId));

        // Telomere extents come from the assembly's chromosomes.tsv, so take the boundary from the
        // reference rather than a constant: contig 0 is the first chromosome of the male genome.
        string firstChrom = _refGen.AllChrNames.First();
        long telomereEnd = _refGen.Telomeres[firstChrom].First().End;

        // Editing immediately inside the boundary leaves the front telomere intact.
        _kar.ApplyInternalDeletion(
            contigId,
            telomereEnd,
            telomereEnd + 1);
        Assert.AreEqual(2, _kar.CountIntactTelomereEnds(contigId));

        // Removing its final base destroys the whole annotation; only the back end remains.
        _kar.ApplyInternalDeletion(
            contigId,
            telomereEnd - 1,
            telomereEnd);
        Assert.AreEqual(1, _kar.CountIntactTelomereEnds(contigId));

        foreach (int id in _kar.ContigIds().Where(id => id != contigId).ToList())
        {
            _kar.ApplyContigDeletion(id);
        }

        var eventParameters = new CNEventPars(CNEventType.TelomereDeletion, 1, 0.01);
        var eventData = Sampling.GenerateCNEventData(_rnd, _kar, eventParameters) as TailEventData;
        Assert.NotNull(eventData);
        Assert.IsFalse(eventData!.Direction);

        long contigLength = _kar.ContigLen(contigId);
        _kar.ApplyInternalInversion(contigId, 0, contigLength);
        eventData = Sampling.GenerateCNEventData(_rnd, _kar, eventParameters) as TailEventData;
        Assert.NotNull(eventData);
        Assert.IsTrue(eventData!.Direction);

        _kar.ApplyTailDeletion(contigId, 1, true);
        Assert.AreEqual(0, _kar.CountIntactTelomereEnds(contigId));
        Assert.IsNull(Sampling.GenerateCNEventData(_rnd, _kar, eventParameters));
    }

    [TestCase(CNEventType.TailDeletion)]
    [TestCase(CNEventType.TailDuplication)]
    [TestCase(CNEventType.TelomereDeletion)]
    [TestCase(CNEventType.TelomereDuplication)]
    [TestCase(CNEventType.InternalDeletion)]
    [TestCase(CNEventType.InternalDuplication)]
    public void TestArmBasedEventRequiresAnInwardCentromere(CNEventType eventType)
    {
        const int contigId = 0;
        foreach (int id in _kar.ContigIds().Where(id => id != contigId).ToList())
        {
            _kar.ApplyContigDeletion(id);
        }
        var centromere = _kar.GetCentromeres(contigId).Single();
        _kar.ApplyInternalDeletion(contigId, centromere.start, centromere.end);
        Assert.AreEqual(2, _kar.CountIntactTelomereEnds(contigId));
        Assert.AreEqual(0, _kar.CountCentromeres(contigId));

        var eventParameters = new CNEventPars(eventType, 1, 0.01);
        Assert.IsNull(Sampling.GenerateCNEventData(_rnd, _kar, eventParameters));
    }

    [TestCase(true)]
    [TestCase(false)]
    public void TestTelomereDeletionStopsAtNearestCentromere(bool direction)
    {
        var ids = _kar.ContigIds().ToList();
        int fusedId = ids[0];
        int donorId = ids[1];
        foreach (int id in ids.Where(id => id != fusedId && id != donorId))
        {
            _kar.ApplyContigDeletion(id);
        }

        FuseContigsForTest(fusedId, donorId);
        Assert.AreEqual(1, _kar.CountContigs());
        Assert.AreEqual(2, _kar.CountCentromeres(fusedId));

        if (direction)
        {
            _kar.ApplyTailDeletion(fusedId, _kar.ContigLen(fusedId) - 1, false);
        }
        else
        {
            _kar.ApplyTailDeletion(fusedId, 1, true);
        }
        Assert.AreEqual(1, _kar.CountIntactTelomereEnds(fusedId));

        var centromeres = _kar.GetCentromeres(fusedId);
        long nearestBoundary = direction
            ? centromeres.Min(centromere => centromere.start)
            : centromeres.Max(centromere => centromere.end);
        var eventParameters = new CNEventPars(CNEventType.TelomereDeletion, 1, 1_000_000);
        var eventData = Sampling.GenerateCNEventData(_rnd, _kar, eventParameters) as TailEventData;
        Assert.NotNull(eventData);
        Assert.AreEqual(direction, eventData!.Direction);
        Assert.AreEqual(nearestBoundary, eventData.Start);

        eventData.ApplyEvent(_kar);
        Assert.AreEqual(2, _kar.CountCentromeres(fusedId));
        Assert.AreEqual(
            direction ? eventData.Length - nearestBoundary : nearestBoundary,
            _kar.ContigLen(fusedId));
    }

    [Test]
    public void TestTailDeletionStopsAtNearestCentromere()
    {
        var ids = _kar.ContigIds().ToList();
        int fusedId = ids[0];
        int donorId = ids[1];
        foreach (int id in ids.Where(id => id != fusedId && id != donorId))
        {
            _kar.ApplyContigDeletion(id);
        }

        FuseContigsForTest(fusedId, donorId);
        Assert.AreEqual(1, _kar.CountContigs());
        var centromeres = _kar.GetCentromeres(fusedId);
        Assert.AreEqual(2, centromeres.Count);

        var observedDirections = new HashSet<bool>();
        var eventParameters = new CNEventPars(CNEventType.TailDeletion, 1, 1_000_000);
        for (int i = 0; i < 100; i++)
        {
            var eventData = Sampling.GenerateCNEventData(_rnd, _kar, eventParameters) as TailEventData;
            Assert.NotNull(eventData);
            observedDirections.Add(eventData!.Direction);
            Assert.AreEqual(
                eventData.Direction
                    ? centromeres.Min(centromere => centromere.start)
                    : centromeres.Max(centromere => centromere.end),
                eventData.Start);
        }
        CollectionAssert.AreEquivalent(new[] { true, false }, observedDirections);
    }

    [Test]
    public void TestTelomereDuplicationUsesTheEligibleContigEnd()
    {
        const int contigId = 0;
        foreach (int id in _kar.ContigIds().Where(id => id != contigId).ToList())
        {
            _kar.ApplyContigDeletion(id);
        }
        _kar.ApplyTailDeletion(contigId, 1, true);
        Assert.AreEqual(1, _kar.CountIntactTelomereEnds(contigId));

        var eventParameters = new CNEventPars(CNEventType.TelomereDuplication, 1, 0.01);
        var eventData = Sampling.GenerateCNEventData(_rnd, _kar, eventParameters) as TailEventData;
        Assert.NotNull(eventData);
        Assert.IsFalse(eventData!.Direction);

        long contigLength = _kar.ContigLen(contigId);
        int contigCount = _kar.CountContigs();
        eventData.ApplyEvent(_kar);
        Assert.AreEqual(contigCount, _kar.CountContigs());
        Assert.AreEqual(
            contigLength + eventData.Length - eventData.Start,
            _kar.ContigLen(contigId));
        Assert.AreEqual(1, _kar.CountIntactTelomereEnds(contigId));
    }
    
    [TestCase(true)]
    [TestCase(false)]
    public void TestBFBEndsWithFinalBreakAndOneCentromere(bool direction)
    {
        long len = _kar.ContigLen(0);
        long removedLength = TEST_FRAC;
        long start = direction ? removedLength : len - removedLength;
        long fusion = len - removedLength;
        long finalBreak = direction ? fusion - TEST_FRAC / 2 : fusion + TEST_FRAC / 2;

        _kar.ApplyBFB(0, start, direction, finalBreak);

        Assert.AreEqual(fusion + TEST_FRAC / 2, _kar.ContigLen(0));
        Assert.AreEqual(1, _kar.CountCentromeres(0));
        Assert.AreEqual(1, _kar.CountIntactTelomereEnds(0));
    }

    [Test]
    public void TestBFBChain()
    {
        foreach (int id in _kar.ContigIds().Where(id => id != 0).ToList())
        {
            _kar.ApplyContigDeletion(id);
        }
        var eventP = new CNEventPars(CNEventType.BreakageFusionBridge, 1, 0.01);
        for (int i = 0; i < 4; i++)
        {
            var eventData = Sampling.GenerateCNEventData(_rnd, _kar, eventP) as BFBEventData;
            Assert.NotNull(eventData);
            Assert.AreEqual(0, eventData!.ContigId);
            eventData.ApplyEvent(_kar);
            Assert.AreEqual(1, _kar.CountCentromeres(0));
        }
    }
    
    [Test]
    public void TestTranslocation()
    {
        long contigLen = _kar.ContigLen(0);
        long chrLen = RegionOps.GetLength(_kar.FindChrRegions("chr1").ToList());
        
        _kar.ApplyTranslocation(0, 1, TEST_FRAC, 2 * TEST_FRAC, true);
        Assert.AreEqual(contigLen + TEST_FRAC, _kar.ContigLen(1));
        Assert.AreEqual(chrLen, RegionOps.GetLength(_kar.FindChrRegions("chr1").ToList()));
        Assert.AreEqual(1, _kar.CountCentromeres(0));
        Assert.AreEqual(1, _kar.CountCentromeres(1));

        _kar.ApplyTranslocation(0, 1, 4 * TEST_FRAC,  3 * TEST_FRAC, true);
        Assert.AreEqual(chrLen, RegionOps.GetLength(_kar.FindChrRegions("chr1").ToList()));
        Assert.AreEqual(contigLen + TEST_FRAC * 2, _kar.ContigLen(0));
        Assert.AreEqual(1, _kar.CountCentromeres(0));
        Assert.AreEqual(1, _kar.CountCentromeres(1));
        Console.WriteLine(_kar);
    }

    [Test]
    public void TestGeneratedTranslocationsKeepBothDerivativesMonocentric()
    {
        var eventP = new CNEventPars(CNEventType.Translocation, 1, 0.1);
        for (int i = 0; i < 200; i++)
        {
            var eventData = Sampling.GenerateCNEventData(_rnd, _kar, eventP) as PairEventData;
            Assert.NotNull(eventData);
            eventData!.ApplyEvent(_kar);
            Assert.AreEqual(1, _kar.CountCentromeres(eventData.ContigIdA));
            Assert.AreEqual(1, _kar.CountCentromeres(eventData.ContigIdB));
        }
    }

    [Test]
    public void TestInvalidTranslocationIsRejectedAtomically()
    {
        var centromereA = _kar.GetCentromeres(0).Single();
        var centromereB = _kar.GetCentromeres(1).Single();
        string before = _kar.ToString();

        Assert.Throws<ArgumentException>(() =>
            _kar.ApplyTranslocation(0, 1, centromereA.start, centromereB.end, false));
        Assert.AreEqual(before, _kar.ToString());
        Assert.AreEqual(1, _kar.CountCentromeres(0));
        Assert.AreEqual(1, _kar.CountCentromeres(1));
    }

    [Test]
    public void TestApplyCNEvent([Values] CNEventType eventType)
    {
        var eventP = new CNEventPars(eventType, 1, 1_000_000, 10);
        Assert.DoesNotThrow(() => ApplyRandomEvent(_rnd, _kar, eventP));
    }

    [Test]
    public void TestRandomEvent([Values] CNEventType eventType, [Values] IntEdgeCases seed)
    {
        var eventP = new CNEventPars(eventType, 1, 1_000_000, 10);
        var eventData = Sampling.GenerateCNEventData(new Random((int) seed), _kar, eventP);
        Assert.NotNull(eventData);
        Assert.DoesNotThrow(() => eventData?.ApplyEvent(_kar));
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
        Assert.AreEqual(_refGen.SexGenomeLen[(int) _kar.Sex], _kar.GenomeLen());
    }
    
    [Test]
    public void TestClean()
    {
        for (int i = 0; i < _refGen.SexGenome[(int) SexType.Female].Count; i++)
        {
            ApplyRandomEvent(_rnd, _kar, _del);
        }
        Assert.AreEqual("[]", _kar.ToString());
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

    [Test]
    public void TestSNV()
    {
        const long loc = 100;
        const int contigID = 0;
        var newNucleotide = Nucleotide.C;

        _kar.ApplyPointMutation(contigID, loc, newNucleotide);
        Assert.AreEqual(46, _kar.CountContigs());
        
        var contig = _kar.GetContig(contigID);
        Assert.AreEqual(1, contig.CountRegions());
        var SNVs = contig.SNVs;
        Assert.NotNull(SNVs);
        Assert.AreEqual(1, SNVs.Count);
        Assert.AreEqual(loc, SNVs[0].Pos);
        Assert.AreEqual(newNucleotide, SNVs[0].Alt);

        // Try a repeated SNV
        newNucleotide = Nucleotide.G;
        _kar.ApplyPointMutation(contigID, loc, newNucleotide);

        Assert.AreEqual(46, _kar.CountContigs());
        
        contig = _kar.GetContig(contigID);
        Assert.AreEqual(1, contig.CountRegions());
        SNVs = contig.SNVs;
        Assert.NotNull(SNVs);
        Assert.AreEqual(1, SNVs.Count);
        Assert.AreEqual(loc, SNVs[0].Pos);
        Assert.AreEqual(newNucleotide, SNVs[0].Alt);
    }

    [Test]
    public void TestSNVWithDeletion()
    {
        const long loc = 100;
        const int contigID = 0;
        var newNucleotide = Nucleotide.C;
        // Apply the SNV
        _kar.ApplyPointMutation(contigID, loc, newNucleotide);
        // Apply a deletion that covers the SNV
        _kar.ApplyInternalDeletion(contigID, 50, 200);
        // Check that the SNV is not present
        var contig = _kar.GetContig(contigID);
        Assert.AreEqual(2, contig.CountRegions());
        var snvs = contig.SNVs;
        Assert.IsEmpty(snvs);
    }

    [Test]
    public void TestSNVWithDuplication()
    {
        const long loc = 100;
        const int contigID = 0;
        var newNucleotide = Nucleotide.C;
        // Apply the SNV
        _kar.ApplyPointMutation(contigID, loc, newNucleotide);
        // Apply a duplication that covers the SNV
        _kar.ApplyInternalDuplication(contigID, 50, 200);
        // Check that the SNV is present in both copies
        var contig = _kar.GetContig(contigID);
        Assert.AreEqual(2, contig.CountRegions());
        var SNVs = contig.SNVs;
        Assert.IsNotEmpty(SNVs);
        Assert.AreEqual(SNVs[0], SNVs[1]);
        Assert.AreEqual(newNucleotide, SNVs[0].Alt);

        // If we alter one region, the other should not be affected
        var secondNucleotide = Nucleotide.G;
        _kar.ApplyPointMutation(contigID, loc, secondNucleotide);
        SNVs = _kar.GetContig(contigID).SNVs;
        Assert.IsNotEmpty(SNVs);
        Assert.AreNotEqual(SNVs[0], SNVs[1]);
        Assert.AreEqual(secondNucleotide, SNVs[0].Alt);
        Assert.AreEqual(newNucleotide, SNVs[1].Alt);
    }

    [Test]
    public void TestCalcChrCopyNumbers()
    {
        var dupEv = new CNEventPars(CNEventType.InternalDuplication, 1, .1);
        int dupCount = 100;
        for (int i = 0; i < dupCount; i++)
        {
            ApplyRandomEvent(_rnd, _kar, dupEv);
        }

        var breaks = _kar.CalcBreaks();
        var chrCopyNumbers = _kar.CalcCNs(breaks);
        // each event should create two new regions
        Assert.AreEqual(24, chrCopyNumbers.Count - dupCount * 2);
    }

    [Test]
    public void TestGetPresentGeneCounts([Values] GeneLT geneType)
    {
        // Assumes _kar is male
        for (int i = 0; i < 23; i++)
        {
            _kar.ApplyContigDeletion(i);
        }
        foreach (var gene in _refGen.SexGeneLists[(int)_kar.Sex][(int) geneType])
        {
            int count = _kar.GeneCounts[(int)geneType][gene.GeneId];
            Assert.AreEqual(gene.Chrom != "chrX" ? 1 : 0, count);
        }
    }
    
    [Test]
    public void TestContigIds()
    {
        for (int i = 0; i < 4; i++)
        {
            _kar.ApplyContigDeletion(i);
        }
        var idsAfter = _kar.ContigIds().ToList();
        Assert.AreEqual(4, idsAfter[0]);
    }

    private double SelectionFraction(CNEventPars ev, int targetContigId, int trials)
    {
        int count = 0;
        for (int i = 0; i < trials; i++)
        {
            var data = Sampling.GenerateCNEventData(_rnd, _kar, ev) as ContigEventData;
            Assert.NotNull(data);
            if (data!.ContigId == targetContigId)
            {
                count++;
            }
        }
        return count / (double) trials;
    }

    // Keeps only the longest and shortest contigs (both single-centromere chromosomes),
    // making length and centromere-count weighting produce visibly different distributions.
    private (int bigId, int smallId) KeepExtremeContigs()
    {
        var ids = _kar.ContigIds().ToList();
        int bigId = ids.OrderByDescending(_kar.ContigLen).First();
        int smallId = ids.OrderBy(_kar.ContigLen).First();
        foreach (int i in ids.Where(i => i != bigId && i != smallId))
        {
            _kar.ApplyContigDeletion(i);
        }
        return (bigId, smallId);
    }

    private double TerminalArmWeight(int contigId, bool front = true, bool back = true)
    {
        var centromeres = _kar.GetCentromeres(contigId);
        long contigLength = _kar.ContigLen(contigId);
        long frontLength = front
            ? centromeres.Select(centromere => centromere.start)
                .Where(length => length > 0)
                .DefaultIfEmpty(0)
                .Min()
            : 0;
        long backLength = back
            ? centromeres.Select(centromere => contigLength - centromere.end)
                .Where(length => length > 0)
                .DefaultIfEmpty(0)
                .Min()
            : 0;
        return frontLength + backLength;
    }

    [Test]
    public void TestInternalSelectionByTerminalArmLength()
    {
        var (bigId, smallId) = KeepExtremeContigs();
        double bigWeight = TerminalArmWeight(bigId);
        double expected = bigWeight / (bigWeight + TerminalArmWeight(smallId));

        var ev = new CNEventPars(CNEventType.InternalDeletion, 1, 0.01);
        double frac = SelectionFraction(ev, bigId, 20000);

        // Picking a terminal arm first weights each contig by the total usable length of its arms.
        Assert.AreEqual(expected, frac, 0.03);
    }

    [TestCase(CNEventType.TailDeletion)]
    [TestCase(CNEventType.TailDuplication)]
    public void TestTailSelectionByTerminalArmLength(CNEventType eventType)
    {
        var (bigId, smallId) = KeepExtremeContigs();
        double bigWeight = TerminalArmWeight(bigId);
        double expected = bigWeight / (bigWeight + TerminalArmWeight(smallId));

        var ev = new CNEventPars(eventType, 1, 0.01);
        Assert.AreEqual(expected, SelectionFraction(ev, bigId, 20000), 0.03);
    }

    [Test]
    public void TestInternalEventsSelectArmsByLengthAndDoNotCrossCentromere()
    {
        const int contigId = 0;
        foreach (int id in _kar.ContigIds().Where(id => id != contigId).ToList())
        {
            _kar.ApplyContigDeletion(id);
        }

        var centromere = _kar.GetCentromeres(contigId).Single();
        long contigLength = _kar.ContigLen(contigId);
        double expectedFront = centromere.start /
            (double) (centromere.start + contigLength - centromere.end);
        int frontCount = 0;
        const int trials = 20000;
        var ev = new CNEventPars(CNEventType.InternalDeletion, 1, 0.1);

        for (int i = 0; i < trials; i++)
        {
            var data = Sampling.GenerateCNEventData(_rnd, _kar, ev) as InternalEventData;
            Assert.NotNull(data);
            if (data!.End <= centromere.start)
            {
                frontCount++;
            }
            else
            {
                Assert.GreaterOrEqual(data.Start, centromere.end);
            }
        }

        Assert.AreEqual(expectedFront, frontCount / (double) trials, 0.03);
    }

    [Test]
    public void TestArmSelectionIgnoresLength()
    {
        var (bigId, _) = KeepExtremeContigs();

        var ev = new CNEventPars(CNEventType.ArmDeletion, 1, 0.1);
        double frac = SelectionFraction(ev, bigId, 20000);

        // Both contigs carry one centromere, so arm events should pick them ~50/50
        // regardless of the large length difference.
        Assert.AreEqual(0.5, frac, 0.03);
    }

    [Test]
    public void TestArmSelectionByCentromereCount()
    {
        // Keep three contigs; fuse two of them into one so it carries two centromeres.
        var ids = _kar.ContigIds().ToList();
        int fusedId = ids[0];
        int donorId = ids[1];
        int singleId = ids[2];
        foreach (int i in ids.Where(i => i != fusedId && i != donorId && i != singleId))
        {
            _kar.ApplyContigDeletion(i);
        }

        // Fuse the donor onto fusedId; the donor is emptied and drops out of the active set.
        FuseContigsForTest(fusedId, donorId);
        Assert.AreEqual(2, _kar.GetCentromeres(fusedId).Count);
        Assert.AreEqual(1, _kar.GetCentromeres(singleId).Count);
        CollectionAssert.AreEquivalent(new[] { fusedId, singleId }, _kar.ContigIds().ToList());

        var ev = new CNEventPars(CNEventType.ArmDeletion, 1, 0.1);
        double frac = SelectionFraction(ev, fusedId, 30000);

        // The fused contig has two centromeres vs one, so it should be picked ~2/3 of the time.
        Assert.AreEqual(2.0 / 3.0, frac, 0.03);
    }

    [TestCase(CNEventType.TelomereDeletion)]
    [TestCase(CNEventType.TelomereDuplication)]
    public void TestTelomereSelectionByEligibleArmLength(CNEventType eventType)
    {
        var (twoEndedId, oneEndedId) = KeepExtremeContigs();
        _kar.ApplyTailDeletion(oneEndedId, 1, true);
        Assert.AreEqual(2, _kar.CountIntactTelomereEnds(twoEndedId));
        Assert.AreEqual(1, _kar.CountIntactTelomereEnds(oneEndedId));

        var ev = new CNEventPars(eventType, 1, 0.01);
        double fraction = SelectionFraction(ev, twoEndedId, 20000);
        double twoEndedWeight = TerminalArmWeight(twoEndedId);
        // The front telomere of oneEndedId was removed, leaving only its back arm eligible.
        double expected = twoEndedWeight /
            (twoEndedWeight + TerminalArmWeight(oneEndedId, front: false));

        Assert.AreEqual(expected, fraction, 0.03);
    }
}
