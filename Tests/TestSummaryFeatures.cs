// Created by Cody Duncan, 2023, codybstrange93@gmail.com

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
public class TestSummaryFeatures
{
    private GenRef _genRef;
    private Karyotype _kar;
    private Random _rnd;
    List<string> _chrs;
    Dictionary<string, List<CopyNumber>> _cnps;
    
    [SetUp]
    public void Setup()
    {
        _genRef = FileIO.GetGenRef(TestIO.HG_19_PATH);
        _rnd = new Random(0);
        _chrs = new List<string> {_genRef.AllChrs[0], _genRef.AllChrs[1]};
        _kar = new Karyotype(_genRef, true);
    }
    [Test]
    public void TestDefaultSegLengths()
    {
        var karXX = new Karyotype(_genRef, true);
        var cn = CopyNumbers.CalcCopyNumbers(_genRef, karXX, true).ToList();
        var cnps = new Dictionary<string, List<CopyNumber>> {{"sample_1", cn}};
        var segLengths = SummaryFeatures.GetSegLengths(cnps);
        Assert.AreEqual(0, segLengths.segs.Count);
        // Count the CN-normal segments of autosomes
        segLengths = SummaryFeatures.GetSegLengths(cnps, true);
        Assert.AreEqual(22, segLengths.segs.Count);
        // Count the CN-normal segments and LoH segments of autosomes
        segLengths = SummaryFeatures.GetSegLengths(cnps, true, true);
        Assert.AreEqual(22, segLengths.segs.Count);
        // Count the CN-normal segments & LoH segments of all chromosomes
        segLengths = SummaryFeatures.GetSegLengths(cnps, true, true, true);
        Assert.AreEqual(23, segLengths.segs.Count);
    }
    [Test]
    public void TestSegLengths()
    {
        var karXX = new Karyotype(_genRef, true);
        // Delete 1000-2000 on chr1
        karXX.ApplyInternalDeletion(0, 1000, 2000);
        // Gain 0-5000 on chr2
        karXX.ApplyInternalDuplication(1, 0, 5000);
        var cn = CopyNumbers.CalcCopyNumbers(_genRef, karXX, true).ToList();
        var cnps = new Dictionary<string, List<CopyNumber>> {{"sample_1", cn}};
        var (segs, max) = SummaryFeatures.GetSegLengths(cnps);
        Assert.AreEqual(2, segs.Count);
        Assert.AreEqual(1000, segs[0]);
        Assert.AreEqual(5000, segs[1]);
        Assert.AreEqual(5000, max);
    }
    [Test]
    public void TestChrCNMatrix()
    {
        var karXX = new Karyotype(_genRef, true);
        var karXY = new Karyotype(_genRef, true);
        karXX.ApplyWGD();
        var cn = CopyNumbers.CalcCopyNumbers(_genRef, karXX, true).ToList();
        var cnps = new Dictionary<string, List<CopyNumber>> {{"sample_1", cn}};
        for (int i = 0; i < 23; i++)
        {
            karXY.ApplyContigDeletion(i);
        }
        cn = CopyNumbers.CalcCopyNumbers(_genRef, karXY, false).ToList();
        cnps.Add("sample_2", cn);
        var autosomes = _genRef.ChrIDsForAutosomes().ToList();
        var matrix = SummaryFeatures.GetChrCopyNumberMatrix(autosomes, cnps);
        // Outer index is chromosome ID
        Assert.AreEqual(22, matrix.Count);
        foreach (var chr in autosomes)
        {
            Assert.AreEqual(2, matrix[chr].Count);
            // WGD karyotype
            Assert.AreEqual(4, matrix[chr]["sample_1"]);
            // Whole-genome halving karyotype
            Assert.AreEqual(1, matrix[chr]["sample_2"]);
        }
    }
    [Test]
    public void TestMKV()
    {
        var karXX = new Karyotype(_genRef, true);
        var karXY = new Karyotype(_genRef, true);
        karXX.ApplyWGD();
        var cn = CopyNumbers.CalcCopyNumbers(_genRef, karXX, true).ToList();
        var cnps = new Dictionary<string, List<CopyNumber>> {{"sample_1", cn}};
        for (int i = 0; i < 23; i++)
        {
            karXY.ApplyContigDeletion(i);
        }
        cn = CopyNumbers.CalcCopyNumbers(_genRef, karXY, false).ToList();
        cnps.Add("sample_2", cn);
        var autosomes = _genRef.ChrIDsForAutosomes().ToList();
        var matrix = SummaryFeatures.GetChrCopyNumberMatrix(autosomes, cnps);
        var mkv = SummaryFeatures.GetMKV(matrix);
        Assert.AreEqual(1.8, mkv, double.Epsilon);
    }
    [Test]
    public void TestAverageAneuploidy()
    {
        var karXX = new Karyotype(_genRef, true);
        var karXY = new Karyotype(_genRef, true);
        karXX.ApplyWGD();
        var cn = CopyNumbers.CalcCopyNumbers(_genRef, karXX, true).ToList();
        var cnps = new Dictionary<string, List<CopyNumber>> {{"sample_1", cn}};
        for (int i = 0; i < 23; i++)
        {
            karXY.ApplyContigDeletion(i);
        }
        cn = CopyNumbers.CalcCopyNumbers(_genRef, karXY, false).ToList();
        cnps.Add("sample_2", cn);
        var autosomes = _genRef.ChrIDsForAutosomes().ToList();
        var matrix = SummaryFeatures.GetChrCopyNumberMatrix(autosomes, cnps);
        var aneuploidy = SummaryFeatures.GetAverageAneuploidy(matrix);
        Assert.AreEqual(0, aneuploidy, double.Epsilon);
    }
    /*
    [Test]
    public void TestDefaultChangepoints()
    {
        var karXX = new Karyotype(_genRef, true);
        var changepoints = SummaryFeatures.GetChangepoints(_genRef, new List<Karyotype> {karXX});
        Assert.AreEqual(0, changepoints.Count);
        // Count the CN-normal segments
        changepoints = SummaryFeatures.GetChangepoints(_genRef, new List<Karyotype> {karXX}, true);
        // Only autosomes are counted
        Assert.AreEqual(22, changepoints.Count);
    }

    [Test]
    public void TestChangepoints()
    {
        var karA = new Karyotype(_genRef, true);
        var karB = new Karyotype(_genRef, true);
        karA.ApplyInternalDeletion(0, 1000, 2000);
        karB.ApplyInternalDuplication(0, 2000, 3000);
        var segs = SummaryFeatures.GetChangepoints(_genRef, new List<Karyotype> {karA, karB});
        Assert.AreEqual(2, segs.Count);
        // Copy-neutral LOH segments are not counted
        karA.ApplyInternalDuplication(23, 1000, 2000);
        segs = SummaryFeatures.GetChangepoints(_genRef, new List<Karyotype> {karA, karB});
        Assert.AreEqual(1, segs.Count);
    }

    [Test]
    public void TestDefaultBreakpointsPerChromosome()
    {
        var karXX = new Karyotype(_genRef, true);
        var breakpoints = SummaryFeatures.GetBreakpointsPerChromosome(_genRef, new List<Karyotype> {karXX});
        Assert.AreEqual(0, breakpoints.Count);
    }

    [Test]
    public void TestBreakpointsPerChromosome()
    {
        var karA = new Karyotype(_genRef, true);
        var karB = new Karyotype(_genRef, true);
        karA.ApplyInternalDeletion(0, 1000, 2000);
        karB.ApplyInternalDuplication(0, 2000, 3000);
        var segs = SummaryFeatures.GetBreakpointsPerChromosome(_genRef, new List<Karyotype> {karA, karB});
        // Two chromosomes have breakpoints: karA H1, and karB H1
        Assert.AreEqual(2, segs.Count);
        // Each chromosome has two breakpoints
        Assert.AreEqual(2, segs[0]);
        Assert.AreEqual(2, segs[1]);
    }

    [Test]
    public void TestGetMinMajCNs()
    {
        var kar = new Karyotype(_genRef, true);
        kar.ApplyInternalDeletion(0, 1000, 2000);
        var cnList = CopyNumbers.CalcCopyNumbers(_genRef, kar, true).ToList();
        var minCNs = SummaryFeatures.GetMinMajCNs(cnList, false);
        // First chromosome
        Assert.AreEqual(1, minCNs[0]);
        Assert.AreEqual(0, minCNs[1]);
        Assert.AreEqual(1, minCNs[2]);
        var majCNs = SummaryFeatures.GetMinMajCNs(cnList, true);
        Assert.AreEqual(1, majCNs[0]);
        Assert.AreEqual(1, majCNs[1]);
        Assert.AreEqual(1, majCNs[2]);
    }
    */

}