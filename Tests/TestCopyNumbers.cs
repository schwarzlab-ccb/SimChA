// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SimChA.DataTypes;
using SimChA.Simulation;
using SimChA.Computation;
using SimChA.IO;

namespace Tests;

[TestFixture]
public class TestCopyNumbers
{
    private Karyotype karyotype;
    private Karyotype referenceMale;
    private Karyotype referenceFemale;
    private Random _rnd;
    private List<CopyNumber> copyNumbers;
    private double ploidy;


    [SetUp]
    public void Setup()
    {
        // karyotype = new Karyotype(false);
        referenceMale = new Karyotype(false);
        referenceFemale = new Karyotype(true);
        _rnd = new Random(0);

        // TODO: modify karyotype
    }

    [Test]
    public void TestCalcPloidyReference()
    {
        var copyNumbersReferenceFemale = CopyNumbers.CalcCopyNumbers(referenceFemale, true).ToList();
        double ploidyReferenceFemale = CopyNumbers.CalcPloidy(copyNumbersReferenceFemale, true);
        Assert.AreEqual(2, ploidyReferenceFemale);

        var copyNumbersReferenceMale = CopyNumbers.CalcCopyNumbers(referenceMale, false).ToList();
        double ploidyReferenceMale = CopyNumbers.CalcPloidy(copyNumbersReferenceMale, false);
        Assert.AreEqual(2, ploidyReferenceMale);
    }

    [Test]
    public void TestCalcPloidy()
    {
        // WGD
        karyotype = new Karyotype(referenceFemale);
        karyotype.ApplyAberration(_rnd, new CNEventP(CNEventType.WholeGenomeDoubling, 1));
        copyNumbers = CopyNumbers.CalcCopyNumbers(karyotype, true).ToList();
        ploidy = CopyNumbers.CalcPloidy(copyNumbers, true);
        Assert.AreEqual(4, ploidy);

        // add a bunch of translocations and inversions and check that ploidy is still 2
        karyotype = new Karyotype(referenceFemale);
        for (int i = 0; i < 100; i++)
        {
            karyotype.ApplyAberration(_rnd, new CNEventP(CNEventType.Translocation, 1));
            var invParams = new Dictionary<string, double> { {"Mean", 0.01} };
            karyotype.ApplyAberration(_rnd, new CNEventP(CNEventType.InternalInversion, 1, invParams));
        }
        copyNumbers = CopyNumbers.CalcCopyNumbers(karyotype, true).ToList();
        ploidy = CopyNumbers.CalcPloidy(copyNumbers, true);
        Assert.AreEqual(2, ploidy);

        // TODO Gain / Loss specific number of chromosomes
    }
}