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
    private float ploidy;


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
        float ploidyReferenceFemale = CopyNumbers.CalcPloidy(copyNumbersReferenceFemale, true);
        Assert.AreEqual(2, ploidyReferenceFemale);

        var copyNumbersReferenceMale = CopyNumbers.CalcCopyNumbers(referenceMale, false).ToList();
        float ploidyReferenceMale = CopyNumbers.CalcPloidy(copyNumbersReferenceMale, false);
        Assert.AreEqual(2, ploidyReferenceMale);
    }

    [Test]
    public void TestCalcPloidy()
    {
        // WGD
        karyotype = new Karyotype(referenceFemale);
        karyotype.ApplyAberration(_rnd, AberrationEnum.WholeGenomeDoubling, new BaseAbbP(1));
        copyNumbers = CopyNumbers.CalcCopyNumbers(karyotype, true).ToList();
        ploidy = CopyNumbers.CalcPloidy(copyNumbers, true);
        Assert.AreEqual(4, ploidy);

        // add a bunch of translocations and inversions and check that ploidy is still 2
        karyotype = new Karyotype(referenceFemale);
        for (int i = 0; i < 100; i++)
        {
            karyotype.ApplyAberration(_rnd, AberrationEnum.Translocation, new BaseAbbP(1));
            karyotype.ApplyAberration(_rnd, AberrationEnum.InternalInversion, new FractionAbbP(1, .01));
        }
        copyNumbers = CopyNumbers.CalcCopyNumbers(karyotype, true).ToList();
        ploidy = CopyNumbers.CalcPloidy(copyNumbers, true);
        Assert.AreEqual(2, ploidy);

        // TODO Gain / Loss specific number of chromosomes
    }


    // [Test]
    // public void TestSplit()
    // {
    //     int remainderLen = _chr1.Length() - 1000;
    //     var rest = _chr1.Split(1000, true);
    //     Assert.AreEqual(1000, _chr1.Length());
    //     Assert.AreEqual(remainderLen, rest.Length());
    // }

    // [Test]
    // public void TestJoin()
    // {
    //     int combinedLen = _chr1.Length() + _chrX.Length();
    //     _chr1.Join(_chrX);
    //     Assert.AreEqual(combinedLen, _chr1.Length());
    // }

    // [Test]
    // public void TestInversion()
    // {
    //     int length = _chr1.Length();
    //     _chr1.InvertRange(length / 4, length * 3 / 4);
    //     Assert.AreEqual(length, _chr1.Length());
    // }

    // [Test]
    // public void TestReplication()
    // {
    //     int length = _chr1.Length() + 900;
    //     _chr1.DuplicateRange(100, 1000);
    //     Assert.AreEqual(length, _chr1.Length());
    // }

    // [Test]
    // public void TestBridge()
    // {
    //     int length = (_chr1.Length() - 1000) * 2;
    //     _chr1.Bridge(1000, true);
    //     Assert.AreEqual(length, _chr1.Length());
    // }

    // [Test]
    // public void TestScatterAndGather()
    // {
    //     int length = _chr1.Length();
    //     _chr1.ScatterAndGather(new List<int> { 1000, 2000, 3000 }, 4, new Random(0));
    //     Assert.AreEqual(length, _chr1.Length());
    // }
}