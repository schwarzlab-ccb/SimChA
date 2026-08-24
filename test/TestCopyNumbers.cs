using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SimChA.Computation;
using SimChA.Data;
using SimChA.EventData;
using SimChA.IO;
using SimChA.Simulation;

namespace Tests;

[TestFixture]
public class TestCopyNumbers
{
    private sealed class SimulatorAccessor : Simulator
    {
        public SimulatorAccessor(Random rnd, RefGen refGen)
            : base(rnd, refGen, new SimParams(), new FitParams(1, 1, 1)) { }

        public static (string Gained, string Lost) CalculateDelta(
            Karyotype before,
            Karyotype after)
            => CalcKaryotypeDiff(before, after);
    }

    private static readonly CNEventType[] DeltaOracleEventTypes =
        Enum.GetValues<CNEventType>()
            .Where(type => type != CNEventType.SNV)
            .ToArray();

    private RefGen _refGen = null!;
    private Random _rnd = null!;
    
    [SetUp]
    public void Setup()
    {
        _refGen = FileIO.ReadGenRef(TestParsing.DATA_PATH, TestParsing.HG_19, TestParsing.GENE_SET);
        _rnd = new Random(0);
    }

    // Preserves the pre-sweep implementation as an independent test oracle.
    private static List<CopyNumber> DiffKaryotypesLegacyOracle(
        Karyotype before,
        Karyotype after)
    {
        var beforeBreaks = before.CalcBreaks();
        var afterBreaks = after.CalcBreaks();
        var jointBreaks = beforeBreaks.Keys
            .Union(afterBreaks.Keys)
            .ToDictionary(
                chrom => chrom,
                chrom => beforeBreaks.GetValueOrDefault(chrom, [])
                    .Concat(afterBreaks.GetValueOrDefault(chrom, []))
                    .ToHashSet()
                    .OrderBy(position => position)
                    .ToList());
        var beforeCNs = CopyNumbers.CalcCNs(before, jointBreaks);
        var afterCNs = CopyNumbers.CalcCNs(after, jointBreaks);

        var deltas = new List<CopyNumber>();
        CopyNumber? run = null;
        for (int i = 0; i < beforeCNs.Count; i++)
        {
            var segment = beforeCNs[i];
            int deltaH1 = afterCNs[i].CNH1 - segment.CNH1;
            int deltaH2 = afterCNs[i].CNH2 - segment.CNH2;
            if (run is not null
                && run.Chrom == segment.Chrom
                && run.End == segment.Start
                && run.CNH1 == deltaH1
                && run.CNH2 == deltaH2)
            {
                run = new CopyNumber(
                    run.Start,
                    segment.End,
                    run.Chrom,
                    deltaH1,
                    deltaH2,
                    0);
            }
            else
            {
                if (run is { CNH1: not 0 } or { CNH2: not 0 })
                {
                    deltas.Add(run);
                }
                run = new CopyNumber(
                    segment.Start,
                    segment.End,
                    segment.Chrom,
                    deltaH1,
                    deltaH2,
                    0);
            }
        }
        if (run is { CNH1: not 0 } or { CNH2: not 0 })
        {
            deltas.Add(run);
        }
        return deltas;
    }

    private static string DeltaSignature(CopyNumber delta)
        => $"{delta.Chrom}:{delta.Start}:{delta.End}:"
           + $"{delta.CNH1}:{delta.CNH2}";

    private static void AssertMatchesLegacyOracle(
        Karyotype before,
        Karyotype after,
        string context)
    {
        var expected = DiffKaryotypesLegacyOracle(before, after)
            .Select(DeltaSignature)
            .ToList();
        var actual = CopyNumbers.DiffKaryotypes(before, after)
            .Select(DeltaSignature)
            .ToList();
        CollectionAssert.AreEqual(expected, actual, context);
    }

    [Test]
    public void TestCalcPloidyReference([Values] SexType sex)
    {
        var kar = new Karyotype(_refGen, sex);
        var cnRef = CopyNumbers.CalcCNs(kar).ToList();
        double ploidyRef = CopyNumbers.CalcPloidy(_refGen, cnRef, sex);
        Assert.AreEqual(2, ploidyRef);
    }

    [Test]
    public void TestCalcPloidyFitness([Values] SexType sex)
    {
        var kar = new Karyotype(_refGen, sex);
        double ploidyRef = SampleStat.CalcPloidy(kar, _refGen);
        Assert.AreEqual(2, ploidyRef);
    }

    [Test]
    public void TestCalcPloidyTetraploid([Values] SexType sex)
    {
        var kar = new Karyotype(_refGen, sex);
        kar.ApplyWGD();
        Assert.AreEqual(sex == SexType.Any ? 88 : 92, kar.CountContigs());
        double tetraploidy = SampleStat.CalcPloidy(kar, _refGen);
        Assert.AreEqual(4, tetraploidy);
    }

    [Test]
    public void TestCalcAutosomeCNs([Values] SexType sex)
    {
        var kar = new Karyotype(_refGen, sex);
        var cnRef = CopyNumbers.CalcCNs(kar).ToList();
        Assert.AreEqual(_refGen.SexChromNames[(int) sex].Count, cnRef.Count);
        double ploidyRef = CopyNumbers.CalcPloidy(_refGen, cnRef, sex);
        Assert.AreEqual(2, ploidyRef);
    }

    [Test]
    public void TestWGSPloidy([Values] SexType sex)
    {
        var kar = new Karyotype(_refGen, sex);
        TestKaryotype.ApplyRandomEvent(_rnd, kar, new CNEventPars(CNEventType.WholeGenomeDoubling, 1));
        var cns = CopyNumbers.CalcCNs(kar).ToList();
        double ploidy = CopyNumbers.CalcPloidy(_refGen, cns, sex);
        Assert.AreEqual(4, ploidy);
        // TODO Gain / Loss specific number of chromosomes
    }

    [Test]
    public void TestNormalPloidy([Values] SexType sex)
    {
        var kar = new Karyotype(_refGen, sex);
        // add a bunch of translocations and inversions and check that ploidy is still 2
        for (int i = 0; i < 3; i++)
        {
            // TestKaryotype.ApplyRandomEvent(_rnd, kar, new CNEventPars(CNEventType.Translocation, 1));
            // TestKaryotype.ApplyRandomEvent(_rnd, kar, new CNEventPars(CNEventType.InternalInversion, 1, 1_000_000));
            kar.ApplyInternalInversion(0, i*10000000, i*20000000);
        }
        var cns = CopyNumbers.CalcCNs(kar).ToList();
        double ploidy = CopyNumbers.CalcPloidy(_refGen, cns, sex);
        Assert.AreEqual(2, ploidy);
    }

    // An internal deletion splits the host region into two, but the copy-number diff must report
    // only the deleted interval as lost, not the split as a loss plus two gains.
    [Test]
    public void TestDiffInternalDeletionReportsOnlyDeletedInterval()
    {
        // Contig 0 is the first autosome on haplotype H1.
        string chrom = _refGen.SexChromNames[(int) SexType.Any][0];
        const long delStart = 1_000_000;
        const long delEnd = 2_000_000;

        var before = new Karyotype(_refGen, SexType.Any);
        var after = new Karyotype(before);
        after.ApplyInternalDeletion(0, delStart, delEnd);

        var deltas = CopyNumbers.DiffKaryotypes(before, after);

        Assert.AreEqual(1, deltas.Count, "exactly one segment should change");
        var delta = deltas[0];
        Assert.AreEqual(chrom, delta.Chrom);
        Assert.AreEqual(delStart, delta.Start);
        Assert.AreEqual(delEnd, delta.End);
        Assert.AreEqual(-1, delta.CNH1, "one copy lost on H1");
        Assert.AreEqual(0, delta.CNH2, "H2 unchanged");
    }

    // A whole-genome doubling gains exactly one copy of every chromosome on both haplotypes,
    // merged into one segment per chromosome.
    [Test]
    public void TestDiffWgdGainsOneCopyGenomeWide()
    {
        var before = new Karyotype(_refGen, SexType.Any);
        var after = new Karyotype(before);
        after.ApplyWGD();

        var deltas = CopyNumbers.DiffKaryotypes(before, after);

        Assert.AreEqual(_refGen.SexChromNames[(int) SexType.Any].Count, deltas.Count);
        Assert.IsTrue(deltas.All(d => d is { CNH1: 1, CNH2: 1 }), "every chromosome gains one copy on both haplotypes");
    }

    [Test]
    public void TestDeltaOutputMergesAdjacentRegionsPerHaplotype()
    {
        string chrom = _refGen.SexChromNames[(int) SexType.Any][0];
        const long start = 1_000_000;
        const long boundary = 2_000_000;
        const long end = 3_000_000;

        var before = new Karyotype(_refGen, SexType.Any);
        var after = new Karyotype(before);
        after.ApplyInternalDuplication(0, start, end);
        after.ApplyInternalDuplication(_refGen.AutosomesCount, boundary, end);

        bool previousPrintDelta = CNEventDesc.PrintDelta;
        try
        {
            CNEventDesc.PrintDelta = true;
            var (gained, lost) = SimulatorAccessor.CalculateDelta(before, after);

            Assert.AreEqual(
                $"[H1:{chrom}[{start}:{end}),H2:{chrom}[{boundary}:{end})]",
                gained);
            Assert.AreEqual("[]", lost);
        }
        finally
        {
            CNEventDesc.PrintDelta = previousPrintDelta;
        }
    }

    [Test]
    public void TestDiffNoChangeIsEmpty()
    {
        var before = new Karyotype(_refGen, SexType.Any);
        var after = new Karyotype(before);
        Assert.IsEmpty(CopyNumbers.DiffKaryotypes(before, after));
    }

    // Inversions preserve copy number, so the copy-number diff reports nothing for them.
    [Test]
    public void TestDiffInversionReportsNothing()
    {
        var before = new Karyotype(_refGen, SexType.Any);
        var after = new Karyotype(before);
        after.ApplyInternalInversion(0, 1_000_000, 2_000_000);
        Assert.IsEmpty(CopyNumbers.DiffKaryotypes(before, after));
    }

    [TestCase(17, SexType.Any)]
    [TestCase(41, SexType.Female)]
    [TestCase(73, SexType.Male)]
    public void TestDiffMatchesLegacyOracleForRandomEventSequences(
        int seed,
        SexType sex)
    {
        const int targetEventCount = 40;
        const int maximumAttempts = 200;
        var random = new Random(seed);
        var current = new Karyotype(_refGen, sex);
        int eventCount = 0;
        int wgdCount = 0;

        for (int attempt = 0;
             attempt < maximumAttempts && eventCount < targetEventCount;
             attempt++)
        {
            CNEventType eventType =
                DeltaOracleEventTypes[random.Next(DeltaOracleEventTypes.Length)];
            if (eventType == CNEventType.WholeGenomeDoubling && wgdCount >= 1)
            {
                continue;
            }

            var eventParameters = new CNEventPars(
                eventType,
                1,
                Frac: 0.1,
                Frag: 4);
            var eventData = Sampling.GenerateCNEventData(
                random,
                current,
                eventParameters);
            if (eventData is null)
            {
                continue;
            }

            var after = new Karyotype(current);
            eventData.ApplyEvent(after);
            AssertMatchesLegacyOracle(
                current,
                after,
                $"seed={seed}, sex={sex}, event={eventCount}, "
                + $"type={eventType}");
            current = after;
            eventCount++;
            if (eventType == CNEventType.WholeGenomeDoubling)
            {
                wgdCount++;
            }
        }

        Assert.AreEqual(
            targetEventCount,
            eventCount,
            $"Could not generate enough events for seed={seed}, sex={sex}");
    }
}
