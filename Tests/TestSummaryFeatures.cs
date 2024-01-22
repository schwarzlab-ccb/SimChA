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
    private Dictionary<string, List<CopyNumber>> GetCNPs(List<Karyotype> kars)
    {
        var cnps = new Dictionary<string, List<CopyNumber>>();
        for(int i = 0; i < kars.Count; i++)
        {
            var kar = kars[i];
            var cn = CopyNumbers.CalcCopyNumbers(_genRef, kar, kar.SexXX).ToList();
            cnps.Add($"sample_{i+1}", cn);
        }
        return cnps;
    }

    [Test]
    public void TestDefaultSegLengths()
    {
        var cnps = GetCNPs(new List<Karyotype> { new(_genRef, true) });
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
        var cnps = GetCNPs(new List<Karyotype> { karXX });
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
        for (int i = 0; i < 23; i++)
        {
            karXY.ApplyContigDeletion(i);
        }
        var cnps = GetCNPs(new List<Karyotype> { karXX, karXY });
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
        for (int i = 0; i < 23; i++)
        {
            karXY.ApplyContigDeletion(i);
        }
        var cnps = GetCNPs(new List<Karyotype> { karXX, karXY });
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
        for (int i = 0; i < 23; i++)
        {
            karXY.ApplyContigDeletion(i);
        }
        var cnps = GetCNPs(new List<Karyotype> { karXX, karXY });
        var autosomes = _genRef.ChrIDsForAutosomes().ToList();
        var matrix = SummaryFeatures.GetChrCopyNumberMatrix(autosomes, cnps);
        var aneuploidy = SummaryFeatures.GetAverageAneuploidy(matrix);
        Assert.AreEqual(0, aneuploidy, double.Epsilon);
    }

    [Test]
    public void TestDefaultChangepoints()
    {
        var karXX = new Karyotype(_genRef, true);
        var cnps = GetCNPs(new List<Karyotype> { karXX });
        var (values, max) = SummaryFeatures.GetChangepointInfo(cnps);
        Assert.AreEqual(0, values.Count);
        Assert.AreEqual(0, max);
        // Count the CN-normal segments
        (values, max) = SummaryFeatures.GetChangepointInfo(cnps, true);
        // Only cn-normal autosomes are counted
        Assert.AreEqual(22, values.Count);
        Assert.AreEqual(0, max);
        // Count CN-normal and LoH autosomes
        (values, max) = SummaryFeatures.GetChangepointInfo(cnps, true, true);
        // Only cn-normal autosomes are counted
        Assert.AreEqual(22, values.Count);
        Assert.AreEqual(0, max);
        // Count all chromosomes
        (values, max) = SummaryFeatures.GetChangepointInfo(cnps, true, true, true);
        // Only cn-normal autosomes are counted
        Assert.AreEqual(23, values.Count);
        Assert.AreEqual(0, max);
    }
    
    [Test]
    public void TestChangepoints()
    {
        var karA = new Karyotype(_genRef, true);
        var karB = new Karyotype(_genRef, true);
        karA.ApplyInternalDeletion(0, 1000, 2000);
        karB.ApplyInternalDuplication(1, 2000, 3000);
        karB.ApplyInternalDuplication(1, 2000, 3000);
        var cnps = GetCNPs(new List<Karyotype> { karA, karB });
        var (values, max) = SummaryFeatures.GetChangepointInfo(cnps);
        Assert.AreEqual(2, values.Count);
        // The step-down on karA, chr1
        Assert.AreEqual(1, values[0]);
        // The step-up on karB, chr2
        Assert.AreEqual(2, values[1]);
        Assert.AreEqual(2, max);
        // Copy-neutral LOH segments are not counted
        karA.ApplyInternalDuplication(23, 1000, 2000);
        cnps = GetCNPs(new List<Karyotype> { karA, karB });
        (values, max) = SummaryFeatures.GetChangepointInfo(cnps);
        Assert.AreEqual(1, values.Count);
        Assert.AreEqual(2, values[0]);
        Assert.AreEqual(2, max);
        // Count the copy-neutral LoH segments (but not CN-normal)
        (values, max) = SummaryFeatures.GetChangepointInfo(cnps, false, true);
        Assert.AreEqual(2, values.Count);
        Assert.AreEqual(0, values[0]);
        Assert.AreEqual(2, values[1]);
        Assert.AreEqual(2, max);
    }
    
    [Test]
    public void TestDefaultBreakpointsPerChromosome()
    {
        var karXX = new Karyotype(_genRef, true);
        var cnps = GetCNPs(new List<Karyotype> { karXX });
        var (values, max) = SummaryFeatures.GetBreakpointsPerChromosome(cnps);
        // Check that there is an entry for each autosome
        Assert.AreEqual(22, values.Count);
        // Check that no breakpoints were found for any of the chromosomes
        foreach(var val in values)
        {
            Assert.AreEqual(0, val);
        }
        Assert.AreEqual(0, max);
    }
    /*
    [Test]
    public void TestBreakpointsPerChromosome()
    {
        var karA = new Karyotype(_genRef, true);
        var karB = new Karyotype(_genRef, true);
        karA.ApplyInternalDeletion(0, 1000, 2000);
        karB.ApplyInternalDuplication(0, 2000, 3000);
        karB.ApplyInternalDuplication(0, 5000, 6000);
        var cnps = GetCNPs(new List<Karyotype> { karA, karB });
        var (values, max) = SummaryFeatures.GetBreakpointsPerChromosome(cnps);
        // Since these are 
        Assert.AreEqual(4, values[0]);
        for (int i = 1; i < values.Count; i++)
        {
            Assert.AreEqual(0, values[i]);
        }
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