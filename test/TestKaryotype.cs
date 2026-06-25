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
    
    [SetUp]
    public void Setup()
    {
        _refGen = FileIO.ReadGenRef(TestParsing.DATA_PATH, TestParsing.HG_19, TestParsing.GENE_SET);
        _kar = new Karyotype(_refGen, SexType.Male);
        _rnd = new Random(0);
        _del = new CNEventPars(CNEventType.ChromDeletion, 1);
    }

    private static List<int[]> CloneGeneCounts(Karyotype kar)
        => kar.GeneCounts.Select(counts => (int[])counts.Clone()).ToList();


    private static List<int[]> CountGenesFromContigs(Karyotype kar)
    {
        var counts = kar.GeneCounts.Select(geneCounts => new int[geneCounts.Length]).ToList();
        foreach (int contigId in kar.ContigIds())
        {
            foreach (Gene gene in kar.GetContig(contigId).Genes)
            {
                counts[(int)gene.ListType][gene.GeneId] += 1;
            }
        }
        return counts;
    }

    private static void AssertGeneCountInvariant(Karyotype kar)
    {
        var observed = CountGenesFromContigs(kar);
        Assert.AreEqual(observed.Count, kar.GeneCounts.Count);
        for (int geneType = 0; geneType < observed.Count; geneType++)
        {
            Assert.AreEqual(observed[geneType].Length, kar.GeneCounts[geneType].Length);
            for (int geneId = 0; geneId < observed[geneType].Length; geneId++)
            {
                Assert.AreEqual(observed[geneType][geneId], kar.GeneCounts[geneType][geneId]);
            }
        }
    }

    private static IEnumerable<CNEventType> GeneratedEventTypes()
        => Enum.GetValues<CNEventType>()
            .Where(eventType => eventType is not CNEventType.Pass and not CNEventType.Skip);

    // Test for each AberrationEnum value
    [Test]
    public void TestWGD()
    {
        long initialGenomeLen = _kar.GenomeLen();
        var initialGeneCounts = CloneGeneCounts(_kar);

        _kar.ApplyWGD();
        Assert.AreEqual(92, _kar.CountContigs());
        Assert.AreEqual(initialGenomeLen * 2, _kar.GenomeLen());
        for (int geneType = 0; geneType < initialGeneCounts.Count; geneType++)
        {
            for (int geneId = 0; geneId < initialGeneCounts[geneType].Length; geneId++)
            {
                Assert.AreEqual(initialGeneCounts[geneType][geneId] * 2, _kar.GeneCounts[geneType][geneId]);
            }
        }
        AssertGeneCountInvariant(_kar);
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
        Assert.That(_kar.ContigIds().Select(_kar.ContigLen), Is.All.GreaterThanOrEqualTo(0));
        AssertGeneCountInvariant(_kar);
    }
    
    [Test]
    public void TestInternalDuplication()
    {
        long len = _kar.ContigLen(0);
        _kar.ApplyInternalDuplication(0, TEST_FRAC, 2 * TEST_FRAC);
        Assert.AreEqual(len + TEST_FRAC, _kar.ContigLen(0));
        AssertGeneCountInvariant(_kar);
        _kar.ApplyInternalDuplication(0, TEST_FRAC, 2 * TEST_FRAC);
        Assert.AreEqual(len + 2 * TEST_FRAC, _kar.ContigLen(0));
        AssertGeneCountInvariant(_kar);
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
        _kar.ApplyTailDuplication(0, TEST_FRAC, true);
        Assert.AreEqual(len, _kar.ContigLen(0));
        Assert.AreEqual(TEST_FRAC, _kar.ContigLen(_kar.ContigIds().Last()));
        _kar.ApplyTailDuplication(0, 2*TEST_FRAC, true);
        Assert.AreEqual(len, _kar.ContigLen(0));
        Assert.AreEqual(2*TEST_FRAC, _kar.ContigLen(_kar.ContigIds().Last()));
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
    public void TestBFBInvertedPalindromicStructure()
    {
        long len = _kar.ContigLen(0);

        _kar.ApplyBFB(0, TEST_FRAC, true);

        var chr1Regions = _kar.GetContig(0).FindChrRegions("chr1").ToList();

        Assert.AreEqual((len - TEST_FRAC) * 2, _kar.ContigLen(0));
        Assert.AreEqual(2, chr1Regions.Count);
        Assert.AreEqual(1, chr1Regions.Count(region => region.Forward));
        Assert.AreEqual(1, chr1Regions.Count(region => !region.Forward));
        Assert.AreEqual(chr1Regions[0].Length, chr1Regions[1].Length);
        AssertGeneCountInvariant(_kar);
    }

    [Test]
    public void TestBFBChain()
    {
        for (int i = 0; i < 4; i++)
        {
            long len = _kar.ContigLen(0);
            _kar.ApplyBFB(0, TEST_FRAC, true);
            Assert.AreEqual((len - TEST_FRAC) * 2, _kar.ContigLen(0));
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

        _kar.ApplyTranslocation(0, 1, 4 * TEST_FRAC,  3 * TEST_FRAC, true);
        Assert.AreEqual(chrLen, RegionOps.GetLength(_kar.FindChrRegions("chr1").ToList()));
        Assert.AreEqual(contigLen + TEST_FRAC * 2, _kar.ContigLen(0));
        Console.WriteLine(_kar);
    }
    
    [Test]
    public void TestTranslocationArmsSwappingConservesLength()
    {
        long genomeLen = _kar.GenomeLen();
        long lenA = _kar.ContigLen(0);
        long lenB = _kar.ContigLen(1);
        long posA = TEST_FRAC;
        long posB = 2 * TEST_FRAC;

        _kar.ApplyTranslocation(0, 1, posA, posB, false);

        Assert.AreEqual(genomeLen, _kar.GenomeLen());
        Assert.AreEqual(posA + (lenB - posB), _kar.ContigLen(0));
        Assert.AreEqual(posB + (lenA - posA), _kar.ContigLen(1));
        AssertGeneCountInvariant(_kar);
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
    public void TestGeneAfterEveryEvent([ValueSource(nameof(GeneratedEventTypes))] CNEventType eventType)
    {
        var eventP = new CNEventPars(eventType, 1, 1_000_000, 10);
        var eventData = Sampling.GenerateCNEventData(new Random(0), _kar, eventP);

        Assert.NotNull(eventData);

        eventData!.ApplyEvent(_kar);

        AssertGeneCountInvariant(_kar);
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

    [Test]
    public void TestInternalSelectionByLength()
    {
        var (bigId, smallId) = KeepExtremeContigs();
        double expected = _kar.ContigLen(bigId) / (double) (_kar.ContigLen(bigId) + _kar.ContigLen(smallId));

        var ev = new CNEventPars(CNEventType.InternalDeletion, 1, 0.01);
        double frac = SelectionFraction(ev, bigId, 20000);

        // Internal events should pick contigs in proportion to their length.
        Assert.AreEqual(expected, frac, 0.03);
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
        _kar.ApplyTranslocation(fusedId, donorId, _kar.ContigLen(fusedId), 0, false);
        Assert.AreEqual(2, _kar.GetCentromeres(fusedId).Count);
        Assert.AreEqual(1, _kar.GetCentromeres(singleId).Count);
        CollectionAssert.AreEquivalent(new[] { fusedId, singleId }, _kar.ContigIds().ToList());

        var ev = new CNEventPars(CNEventType.ArmDeletion, 1, 0.1);
        double frac = SelectionFraction(ev, fusedId, 30000);

        // The fused contig has two centromeres vs one, so it should be picked ~2/3 of the time.
        Assert.AreEqual(2.0 / 3.0, frac, 0.03);
    }
}
