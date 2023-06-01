// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SimChA.DataTypes;
using SimChA.Simulation;
using SimChA.Computation;
using SimChA.EventData;

namespace Tests;

[TestFixture]
public class TestCopyNumbers
{
    private Random _rnd;
    
    [SetUp]
    public void Setup()
    {
        _rnd = new Random(0);
        HGRef.Assembly = GenomeAssembly.hg19;
    }

    [Test]
    public void TestCalcPloidyReference([Values] bool sex)
    {
        var kar = new Karyotype(sex);
        var cnRef = CopyNumbers.CalcCopyNumbers(kar, sex).ToList();
        double ploidyRef = CopyNumbers.CalcPloidy(cnRef, sex);
        Assert.AreEqual(2, ploidyRef);
    }

    [Test]
    public void TestWGSPloid([Values] bool sex)
    {
        var kar = new Karyotype(sex);
        TestKaryotype.ApplyRandomEvent(_rnd, kar, new CNEventPars(CNEventType.WholeGenomeDoubling, 1));
        var cns = CopyNumbers.CalcCopyNumbers(kar, sex).ToList();
        var ploidy = CopyNumbers.CalcPloidy(cns, sex);
        Assert.AreEqual(4, ploidy);
        // TODO Gain / Loss specific number of chromosomes
    }

    [Test]
    public void TestNormalPloidy([Values] bool sex)
    {
        var kar = new Karyotype(sex);
        // add a bunch of translocations and inversions and check that ploidy is still 2
        for (int i = 0; i < 100; i++)
        {
            TestKaryotype.ApplyRandomEvent(_rnd, kar, new CNEventPars(CNEventType.Translocation, 1));
            var invParams = new Dictionary<string, double> { {"Mean", 0.01} };
            TestKaryotype.ApplyRandomEvent(_rnd, kar, new CNEventPars(CNEventType.InternalInversion, 1, invParams));
        }
        var cns = CopyNumbers.CalcCopyNumbers(kar, sex).ToList();
        var ploidy = CopyNumbers.CalcPloidy(cns, sex);
        Assert.AreEqual(2, ploidy);
    }

    [Test]
    public void TestDefaultSegPoints()
    {
        var karXX = new Karyotype(true);
        var karXY = new Karyotype(false);
        var segs = CopyNumbers.GetSegPoints(new List<Karyotype> {karXX, karXY});
        foreach (var seg in segs)
        {
            Assert.AreEqual(0, seg.Value.First());
            Assert.AreEqual(HGRef.GetChromLen(seg.Key), seg.Value.Last());
        }   
    }
    
    [Test]
    public void TestCutSegPoints()
    {
        var karXX = new Karyotype(true);
        var karXY = new Karyotype(false);
        karXX.ApplyInternalDeletion(0, 1000, 2000);
        karXY.ApplyInternalDeletion(0, 2000, 3000);
        var segs = CopyNumbers.GetSegPoints(ChrNo.chr1, new List<Karyotype> {karXX, karXY});
        var expected = new List<long> {0, 1000, 2000, 3000, HGRef.GetChromLen(ChrNo.chr1)};
        Assert.AreEqual(expected, segs);
    }
    
    [Test]
    public void TestCutCNs()
    {
        var karXX = new Karyotype(true);
        karXX.ApplyInternalDeletion(0, 1000, 2000);
        var segs = CopyNumbers.GetSegPoints(ChrNo.chr1, new List<Karyotype> {karXX});
        var cns = CopyNumbers.CalcChrCopyNumbers(karXX.FindRegionsOfChr(ChrNo.chr1).ToList(), karXX.GetMissingOfChr(ChrNo.chr1), segs, ChrNo.chr1);
        Console.WriteLine(string.Join(", ", cns));
    }
}