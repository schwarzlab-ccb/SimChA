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
}