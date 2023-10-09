// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SimChA.DataTypes;
using SimChA.Simulation;
using SimChA.Computation;
using SimChA.EventData;
using SimChA.IO;

namespace Tests;

[TestFixture]
public class TestCopyNumbers
{
    private GenRef _genRef;
    private Karyotype _kar;
    private Random _rnd;
    
    [SetUp]
    public void Setup()
    {
        _genRef = FileIO.GetGenRef(TestIO.HG_19_PATH);
        _rnd = new Random(0);
    }

    [Test]
    public void TestCalcPloidyReference([Values] bool sex)
    {
        _kar = new Karyotype(_genRef, sex);
        var cnRef = CopyNumbers.CalcCopyNumbers(_genRef, _kar, sex).ToList();
        double ploidyRef = CopyNumbers.CalcPloidy(_genRef, cnRef, sex);
        Assert.AreEqual(2, ploidyRef);
    }

    [Test]
    public void TestWGSPloid([Values] bool sex)
    {
        _kar = new Karyotype(_genRef, sex);
        TestKaryotype.ApplyRandomEvent(_rnd, _kar, new CNEventPars(CNEventType.WholeGenomeDoubling, 1));
        var cns = CopyNumbers.CalcCopyNumbers(_genRef, _kar, sex).ToList();
        double ploidy = CopyNumbers.CalcPloidy(_genRef, cns, sex);
        Assert.AreEqual(4, ploidy);
        // TODO Gain / Loss specific number of chromosomes
    }

    [Test]
    public void TestNormalPloidy([Values] bool sex)
    {
        _kar = new Karyotype(_genRef, sex);
        // add a bunch of translocations and inversions and check that ploidy is still 2
        for (int i = 0; i < 100; i++)
        {
            TestKaryotype.ApplyRandomEvent(_rnd, _kar, new CNEventPars(CNEventType.Translocation, 1));
            TestKaryotype.ApplyRandomEvent(_rnd, _kar, new CNEventPars(CNEventType.InternalInversion, 1, 1_000_000));
        }
        var cns = CopyNumbers.CalcCopyNumbers(_genRef, _kar, sex).ToList();
        var ploidy = CopyNumbers.CalcPloidy(_genRef, cns, sex);
        Assert.AreEqual(2, ploidy);
    }

    [Test]
    public void TestDefaultSegPoints()
    {
        var karXX = new Karyotype(_genRef, true);
        var karXY = new Karyotype(_genRef, false);
        var segs = CopyNumbers.GetSegPoints(_genRef, new List<Karyotype> {karXX, karXY});
        foreach (var seg in segs)
        {
            Assert.AreEqual(0, seg.Value.First());
            Assert.AreEqual(_genRef.ChrLengths[seg.Key], seg.Value.Last());
        }   
    }

    [Test]
    public void TestDefaultSegLengths()
    {
        var karXX = new Karyotype(_genRef, true);
        var karXY = new Karyotype(_genRef, false);
        var segs = CopyNumbers.GetSegLengths(_genRef, new List<Karyotype> {karXX, karXY});
        foreach (var seg in segs)
        {
            Assert.AreEqual(_genRef.ChrLengths[seg.Key], seg.Value.First());
        }
    }

    [Test]
    public void TestSegLengths()
    {
        var karA = new Karyotype(_genRef, true);
        var karB = new Karyotype(_genRef, true);
        karA.ApplyInternalDeletion(0, 1000, 2000);
        karB.ApplyInternalDuplication(0, 2000, 3000);
        var segs = CopyNumbers.GetSegLengths(_genRef, new List<Karyotype> {karA, karB});
        // The duplicated and deleted region lengths for haplotype A
        Assert.AreEqual(1000, segs["chr1"][0]);
        Assert.AreEqual(1000, segs["chr1"][1]);
        // Undisturbed beginning section of karB, haplotype A
        Assert.AreEqual(2000, segs["chr1"][2]);
        // Rest of haplotype A for karA and karB
        Assert.AreEqual(_genRef.ChrLengths["chr1"]-2000, segs["chr1"][3]);
        Assert.AreEqual(_genRef.ChrLengths["chr1"]-2000, segs["chr1"][4]);
        // Haplotype B for the two different karyotypes
        Assert.AreEqual(_genRef.ChrLengths["chr1"], segs["chr1"][5]);
        Assert.AreEqual(_genRef.ChrLengths["chr1"], segs["chr1"][6]);
    }
    
    [Test]
    public void TestCutSegPoints()
    {
        var karXX = new Karyotype(_genRef, true);
        var karXY = new Karyotype(_genRef, false);
        karXX.ApplyInternalDeletion(0, 1000, 2000);
        karXY.ApplyInternalDeletion(0, 2000, 3000);
        var segs = CopyNumbers.GetSegPoints(_genRef, "chr1", new List<Karyotype> {karXX, karXY});
        var expected = new List<long> {0, 1000, 2000, 3000, _genRef.ChrLengths["chr1"]};
        Assert.AreEqual(expected, segs);
    }
    
    [Test]
    public void TestCutCNs()
    {
        var karXX = new Karyotype(_genRef, true);
        karXX.ApplyInternalDeletion(0, 1000, 2000);
        string chrNo = "chr1";
        var segs = CopyNumbers.GetSegPoints(_genRef, chrNo, new List<Karyotype> {karXX});
        var cns = CopyNumbers.CalcChrCopyNumbers(
            karXX.FindRegionsOfChr(chrNo).ToList(), 
            karXX.GetMissingOfChr(chrNo), 
            segs, 
            chrNo, 
            false);
        Console.WriteLine(string.Join(", ", cns));
    }
}