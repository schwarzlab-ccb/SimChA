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
    
    [SetUp]
    public void Setup()
    {
        _genRef = FileIO.GetGenRef(TestIO.HG_19_PATH);
        _rnd = new Random(0);
    }

    /*[Test]
    public void TestDefaultSegLengths()
    {
        var karXX = new Karyotype(_genRef, true);
        var segLengths = SummaryFeatures.GetSegLengths(_genRef, new List<Karyotype> {karXX});
        Assert.AreEqual(0, segLengths.Count);
        // Count the CN-normal segments
        segLengths = SummaryFeatures.GetSegLengths(_genRef, new List<Karyotype> {karXX}, true);
        // Only autosomes are counted
        Assert.AreEqual(22, segLengths.Count);
    }

    [Test]
    public void TestSegLengths()
    {
        var karA = new Karyotype(_genRef, true);
        var karB = new Karyotype(_genRef, true);
        karA.ApplyInternalDeletion(0, 1000, 2000);
        karB.ApplyInternalDuplication(0, 2000, 3000);
        
        var segs = SummaryFeatures.GetSegLengths(_genRef, new List<Karyotype> {karA, karB});
        Assert.AreEqual(2, segs.Count);
        // Copy-neutral LOH segments are not counted
        karA.ApplyInternalDuplication(23, 1000, 2000);
        segs = SummaryFeatures.GetSegLengths(_genRef, new List<Karyotype> {karA, karB});
        Assert.AreEqual(1, segs.Count);
    }

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