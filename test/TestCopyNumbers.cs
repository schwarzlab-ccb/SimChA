using System;
using System.Linq;
using NUnit.Framework;
using SimChA.Computation;
using SimChA.Data;
using SimChA.EventData;
using SimChA.IO;

namespace Tests;

[TestFixture]
public class TestCopyNumbers
{
    private RefGen _refGen = null!;
    private Random _rnd = null!;
    
    [SetUp]
    public void Setup()
    {
        _refGen = FileIO.ReadGenRef(TestParsing.DATA_PATH, TestParsing.HG_19, TestParsing.GENE_SET);
        _rnd = new Random(0);
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
}