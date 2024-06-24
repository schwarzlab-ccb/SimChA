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
    public void TestDefaultMeanSegLength()
    {
        var cnps = GetCNPs(new List<Karyotype> { new(_genRef, true) });
        var meanSegLengths = SummaryFeatures.GetMeanSegLength(cnps);
        // One genome returned with 0 average segment length
        Assert.AreEqual(1, meanSegLengths.Count);
        Assert.AreEqual(0, meanSegLengths[0]);
        // Count the CN-normal segments of autosomes
        meanSegLengths = SummaryFeatures.GetMeanSegLength(cnps, true);
        Assert.AreEqual(1, meanSegLengths.Count);
        Assert.AreEqual(_genRef.AutosomeLinLen/22.0, meanSegLengths[0]);
        // Count the CN-normal segments and LoH segments of autosomes
        meanSegLengths = SummaryFeatures.GetMeanSegLength(cnps, true, true);
        Assert.AreEqual(1, meanSegLengths.Count);
        Assert.AreEqual(_genRef.AutosomeLinLen/22.0, meanSegLengths[0]);
        // Count the CN-normal segments & LoH segments of all chromosomes
        meanSegLengths = SummaryFeatures.GetMeanSegLength(cnps, true, true, true);
        Assert.AreEqual(1, meanSegLengths.Count);
        Assert.AreEqual(_genRef.GetGenomeLen(true, false)/23.0, meanSegLengths[0]);
        // Mean segment length weighted by copy-number
        meanSegLengths = SummaryFeatures.GetMeanSegLength(cnps, true, true, true, true);
        Assert.AreEqual(1, meanSegLengths.Count);
        Assert.AreEqual(_genRef.GetGenomeLen(true, false)/23.0, meanSegLengths[0]);
    }

    [Test]
    public void TestMeanSegLength()
    {
        var karXX = new Karyotype(_genRef, true);
        // make the genome haploid, except 4 copies of chr1-h1, and 2 copies of chr2-h1
        for (int i = 0; i < 3; i++)
        {
            karXX.ApplyContigDuplication(0);
        }
        karXX.ApplyContigDuplication(1);
        karXX.ApplyContigDeletion(24);
        for (int i = 2; i < 23; i++)
        {
            karXX.ApplyContigDeletion(i);
        }
        var cnps = GetCNPs(new List<Karyotype> { karXX });
        var meanSegLengths = SummaryFeatures.GetMeanSegLength(cnps);
        // One genome returned with 0 average segment length
        Assert.AreEqual(1, meanSegLengths.Count);
        // chr2 should be missing since it is an LOH region 
        // (should also be missing from numerator in average, i.e. 21 instead of 22)
        double expectedLen = _genRef.AutosomeLinLen - _genRef.ChrLengths[_genRef.AllChrs[1]];
        Assert.AreEqual(expectedLen/21.0, meanSegLengths[0]);
        // Count the CN-normal segments of autosomes
        // Same situation as above, i.e. chr2 is missing
        meanSegLengths = SummaryFeatures.GetMeanSegLength(cnps, true);
        Assert.AreEqual(1, meanSegLengths.Count);
        Assert.AreEqual(expectedLen/21.0, meanSegLengths[0], double.Epsilon);
        // Count the CN-normal segments and LoH segments of autosomes
        // Chr2 now counted
        meanSegLengths = SummaryFeatures.GetMeanSegLength(cnps, true, true);
        Assert.AreEqual(1, meanSegLengths.Count);
        Assert.AreEqual(_genRef.AutosomeLinLen/22.0, meanSegLengths[0], double.Epsilon);
        // Count the CN-normal segments & LoH segments of all chromosomes
        meanSegLengths = SummaryFeatures.GetMeanSegLength(cnps, true, true, true);
        Assert.AreEqual(1, meanSegLengths.Count);
        Assert.AreEqual(_genRef.GetGenomeLen(true, false)/23.0, meanSegLengths[0], double.Epsilon);
        // Count the mean segment length weighted by copy-number
        meanSegLengths = SummaryFeatures.GetMeanSegLength(cnps, true, true, true, true);
        Assert.AreEqual(1, meanSegLengths.Count);
        // There are 5 copies of chr1, 2 copies of chr2, and one copy of all other chromosomes (including chrX)
        // i.e. 28 'chromosomes' in total
        expectedLen = _genRef.GetGenomeLen(true, false) + 4.0*_genRef.ChrLengths[_genRef.AllChrs[0]] + _genRef.ChrLengths[_genRef.AllChrs[1]];
        Assert.AreEqual(expectedLen/28.0, meanSegLengths[0], double.Epsilon);
    }

    [Test]
    public void TestDefaultGetBreakpoints()
    {
        var cnps = GetCNPs(new List<Karyotype> { new(_genRef, true) });
        // Let there be one event to stop division by 0
        var eventCounts = new Dictionary<string, int> {{"sample_1", 1}};
        // Autosomes only
        var breakpoints = SummaryFeatures.GetBreakpoints(cnps, eventCounts, false);
        Assert.AreEqual(1, breakpoints.Count);
        Assert.AreEqual(0, breakpoints[0]);
        // Including sex chromosomes
        breakpoints = SummaryFeatures.GetBreakpoints(cnps, eventCounts, true);
        Assert.AreEqual(1, breakpoints.Count);
        Assert.AreEqual(0, breakpoints[0]);
    }

    [Test]
    public void TestGetBreakpoints()
    {
        var karXX = new Karyotype(_genRef, true);
        // Delete 1000-2000 on chr1
        karXX.ApplyInternalDeletion(0, 1000, 2000);
        // Gain 0-5000 on chr2
        karXX.ApplyInternalDuplication(1, 0, 5000);
        // Gain on chrX
        karXX.ApplyInternalDuplication(22, 1000, 2000);
        var cnps = GetCNPs(new List<Karyotype> { karXX });
        // Let there be 3 events
        var eventCounts = new Dictionary<string, int> {{"sample_1", 3}};
        // Autosomes only
        var breakpoints = SummaryFeatures.GetBreakpoints(cnps, eventCounts, false);
        Assert.AreEqual(1, breakpoints.Count);
        Assert.AreEqual(1, breakpoints[0]);
        // Including sex chromosomes
        breakpoints = SummaryFeatures.GetBreakpoints(cnps, eventCounts, true);
        Assert.AreEqual(1, breakpoints.Count);
        Assert.AreEqual(5.0/3.0, breakpoints[0]);
    }

    [Test]
    public void TestDefaultStratifiedSegLengths()
    {
        var cnps = GetCNPs(new List<Karyotype> { new(_genRef, true) });
        var stratifiedSegLengths = SummaryFeatures.GetStratifiedSegLengths(cnps);
        // Should have three entries corresponding to cn < 2, cn == 2, cn > 2
        Assert.AreEqual(3, stratifiedSegLengths.Count);
        
        Assert.AreEqual(0, stratifiedSegLengths[0].segs.Count);
        Assert.AreEqual(22, stratifiedSegLengths[1].segs.Count);
        Assert.AreEqual(0, stratifiedSegLengths[2].segs.Count);
        // Weights are by default equal
        Assert.AreEqual(1.0/3, stratifiedSegLengths[0].weight, double.Epsilon);
        Assert.AreEqual(1.0/3, stratifiedSegLengths[1].weight, double.Epsilon);
        Assert.AreEqual(1.0/3, stratifiedSegLengths[2].weight, double.Epsilon);
        // Weights not equal
        stratifiedSegLengths = SummaryFeatures.GetStratifiedSegLengths(cnps, true);
        Assert.AreEqual(3, stratifiedSegLengths.Count);
        
        Assert.AreEqual(0, stratifiedSegLengths[0].segs.Count);
        Assert.AreEqual(22, stratifiedSegLengths[1].segs.Count);
        Assert.AreEqual(0, stratifiedSegLengths[2].segs.Count);
        // Weights now calculated according to number of segments
        Assert.AreEqual(0, stratifiedSegLengths[0].weight, double.Epsilon);
        Assert.AreEqual(1, stratifiedSegLengths[1].weight, double.Epsilon);
        Assert.AreEqual(0, stratifiedSegLengths[2].weight, double.Epsilon);
    }

    [Test]
    public void TestStratifiedSegLengths()
    {
        var karXX = new Karyotype(_genRef, true);
        // Delete 1000-2000 on chr1
        karXX.ApplyInternalDeletion(0, 1000, 2000);
        // Gain 0-5000 on chr2
        karXX.ApplyInternalDuplication(1, 0, 5000);
        var cnps = GetCNPs(new List<Karyotype> { karXX });
        var stratifiedSegLengths = SummaryFeatures.GetStratifiedSegLengths(cnps);
        Assert.AreEqual(3, stratifiedSegLengths.Count);
        // cn < 2
        Assert.AreEqual(1, stratifiedSegLengths[0].segs.Count);
        Assert.AreEqual(1000, stratifiedSegLengths[0].segs[0]);
        // cn == 2
        Assert.AreEqual(23, stratifiedSegLengths[1].segs.Count);
        // cn == 3
        Assert.AreEqual(1, stratifiedSegLengths[2].segs.Count);
        Assert.AreEqual(5000, stratifiedSegLengths[2].segs[0]);

        // Weights should all be equal by default
        Assert.AreEqual(1.0/3, stratifiedSegLengths[0].weight, double.Epsilon);
        Assert.AreEqual(1.0/3, stratifiedSegLengths[1].weight, double.Epsilon);
        Assert.AreEqual(1.0/3, stratifiedSegLengths[2].weight, double.Epsilon);

        // Recalculating weights based off count
        stratifiedSegLengths = SummaryFeatures.GetStratifiedSegLengths(cnps, true);
        Assert.AreEqual(3, stratifiedSegLengths.Count);
        // cn < 2
        Assert.AreEqual(1, stratifiedSegLengths[0].segs.Count);
        Assert.AreEqual(1000, stratifiedSegLengths[0].segs[0]);
        // cn == 2
        Assert.AreEqual(23, stratifiedSegLengths[1].segs.Count);
        // cn > 2
        Assert.AreEqual(1, stratifiedSegLengths[2].segs.Count);
        Assert.AreEqual(5000, stratifiedSegLengths[2].segs[0]);
        // Weights now calculated according to number of segments
        Assert.AreEqual(1.0/25, stratifiedSegLengths[0].weight, double.Epsilon);
        Assert.AreEqual(23.0/25, stratifiedSegLengths[1].weight, double.Epsilon);
        Assert.AreEqual(1.0/25, stratifiedSegLengths[2].weight, double.Epsilon);
    }

    [Test]
    public void TestDefaultSegLengths()
    {
        var cnps = GetCNPs(new List<Karyotype> { new(_genRef, true) });
        var noCutoff = -1;
        var segLengths = SummaryFeatures.GetSegLengths(cnps, noCutoff);
        Assert.AreEqual(0, segLengths.Count);
        // Count the CN-normal segments of autosomes
        segLengths = SummaryFeatures.GetSegLengths(cnps, noCutoff, true);
        Assert.AreEqual(22, segLengths.Count);
        // Count the CN-normal segments and LoH segments of autosomes
        segLengths = SummaryFeatures.GetSegLengths(cnps, noCutoff, true, true);
        Assert.AreEqual(22, segLengths.Count);
        // Count the CN-normal segments & LoH segments of all chromosomes
        segLengths = SummaryFeatures.GetSegLengths(cnps, noCutoff, true, true, true);
        Assert.AreEqual(23, segLengths.Count);
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
        var segs = SummaryFeatures.GetSegLengths(cnps);
        Assert.AreEqual(2, segs.Count);
        Assert.AreEqual(1000, segs[0]);
        Assert.AreEqual(5000, segs[1]);
        Assert.AreEqual(5000, segs.Max());
    }

    [Test]
    public void TestDefaultSegLengthsCutoff()
    {
        var cnps = GetCNPs(new List<Karyotype> { new(_genRef, true) });
        // Cutoff of set to exclude all chromosomes
        var cutoff = 10;
        var segLengths = SummaryFeatures.GetSegLengths(cnps, cutoff);
        Assert.AreEqual(0, segLengths.Count);
        segLengths = SummaryFeatures.GetSegLengths(cnps, cutoff, true);
        Assert.AreEqual(0, segLengths.Count);
        // Count the CN-normal segments and LoH segments of autosomes
        segLengths = SummaryFeatures.GetSegLengths(cnps, cutoff, true, true);
        Assert.AreEqual(0, segLengths.Count);
        // Count the CN-normal segments & LoH segments of all chromosomes
        segLengths = SummaryFeatures.GetSegLengths(cnps, cutoff, true, true, true);
        Assert.AreEqual(0, segLengths.Count);
    }

    [Test]
    public void TestSegLengthsCutoff()
    {
        var karXX = new Karyotype(_genRef, true);
        // Delete 100_000_000-200_000_000 on chr1
        karXX.ApplyInternalDeletion(0, 100_000_000, 200_000_000);
        var cutoff = 20_000_000;
        var cnps = GetCNPs(new List<Karyotype> { karXX });
        var segs = SummaryFeatures.GetSegLengths(cnps, cutoff);
        Assert.AreEqual(0, segs.Count);
        long noCutoff = -1;
        segs = SummaryFeatures.GetSegLengths(cnps, noCutoff);
        Assert.AreEqual(1, segs.Count);
        Assert.AreEqual(100_000_000, segs[0], double.Epsilon);
    }

    [Test]
    public void TestDefaultGetMeanPloidy()
    {
        var karXX = new Karyotype(_genRef, true);
        var karXY = new Karyotype(_genRef, false);
        var cnps = GetCNPs(new List<Karyotype> { karXX, karXY });
        var sexDict = new Dictionary<string, bool> {{cnps.Keys.ToList()[0], true }, {cnps.Keys.ToList()[1], false}};
        var meanPloidy = SummaryFeatures.GetMeanPloidy(_genRef, cnps, sexDict);
        Assert.AreEqual(2, meanPloidy);
    }

    [Test]
    public void TestGetMeanPloidy()
    {
        var karXX = new Karyotype(_genRef, true);
        // WGD and deletion of two copies of X chromosome
        karXX.ApplyWGD();
        karXX.ApplyContigDeletion(22);
        karXX.ApplyContigDeletion(45);
        var karXY = new Karyotype(_genRef, false);
        // Delete all contigs from chr1-22 (i.e. one copy of each autosome)
        for (int i = 0; i < 22; i++)
        {
            karXY.ApplyContigDeletion(i);
        }
        var cnps = GetCNPs(new List<Karyotype> { karXX, karXY });
        var sexDict = new Dictionary<string, bool> {{cnps.Keys.ToList()[0], true }, {cnps.Keys.ToList()[1], false}};
        // Autosomes only
        var mean = SummaryFeatures.GetMeanPloidy(_genRef, cnps, sexDict);
        Assert.AreEqual(2.5, mean);
    }

    [Test]
    public void TestDefaultGetPloidy()
    {
        var karXX = new Karyotype(_genRef, true);
        var karXY = new Karyotype(_genRef, false);
        var cnps = GetCNPs(new List<Karyotype> { karXX, karXY });
        var sexDict = new Dictionary<string, bool> {{cnps.Keys.ToList()[0], true }, {cnps.Keys.ToList()[1], false}};
        // Autosomes only
        var ploidies = SummaryFeatures.GetPloidy(_genRef, cnps, sexDict);
        Assert.AreEqual(2, ploidies.Count);
        Assert.AreEqual(2, ploidies[0]);
        Assert.AreEqual(2, ploidies[1]);
        // Sex Chromosomes included
        ploidies = SummaryFeatures.GetPloidy(_genRef, cnps, sexDict, true);
        Assert.AreEqual(2, ploidies.Count);
        Assert.AreEqual(2, ploidies[0]);
        Assert.AreEqual(2, ploidies[1]);
    }

    [Test]
    public void TestGetPloidy()
    {
        var karXX = new Karyotype(_genRef, true);
        // WGD and deletion of two copies of X chromosome
        karXX.ApplyWGD();
        karXX.ApplyContigDeletion(22);
        karXX.ApplyContigDeletion(45);
        var karXY = new Karyotype(_genRef, false);
        // Delete all contigs from chr1-22 (i.e. one copy of each autosome)
        for (int i = 0; i < 22; i++)
        {
            karXY.ApplyContigDeletion(i);
        }
        var cnps = GetCNPs(new List<Karyotype> { karXX, karXY });
        var sexDict = new Dictionary<string, bool> {{cnps.Keys.ToList()[0], true }, {cnps.Keys.ToList()[1], false}};
        // Autosomes only
        var ploidies = SummaryFeatures.GetPloidy(_genRef, cnps, sexDict);
        Assert.AreEqual(2, ploidies.Count);
        Assert.AreEqual(4, ploidies[0]);
        Assert.AreEqual(1, ploidies[1]);
        // Sex Chromosomes included
        ploidies = SummaryFeatures.GetPloidy(_genRef, cnps, sexDict, true);
        Assert.AreEqual(2, ploidies.Count);
        var expectedXX = (4.0*_genRef.GetGenomeLen(true, false) - 2.0*_genRef.ChrLengths[_genRef.AllChrs[22]])/_genRef.GetGenomeLen(true, false);
        Assert.AreEqual(expectedXX, ploidies[0], 1e-6);
        var expectedXY = 2.0*(_genRef.GetGenomeLen(false) - _genRef.AutosomeLinLen)/_genRef.GetGenomeLen(false);
        Assert.AreEqual(expectedXY, ploidies[1], 1e-6);
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
    public void TestDefaultMeanBreakpoints()
    {
        var karXX = new Karyotype(_genRef, true);
        var cnps = GetCNPs(new List<Karyotype> { karXX });
        var eventCounts = new Dictionary<string, int> {{cnps.Keys.ToList()[0], 1 }};
        var mean = SummaryFeatures.GetMeanBreakpoints(cnps, eventCounts, false);
        Assert.AreEqual(0, mean);
    }

    [Test]
    public void TestMeanBreakpoints()
    {
        var karXX = new Karyotype(_genRef, true);
        // Missing segment from the start
        karXX.ApplyInternalDeletion(0, 0, 1000);
        karXX.ApplyInternalDuplication(1, 1000, 2000);
        var cnps = GetCNPs(new List<Karyotype> { karXX });
        var nEvents = 2;
        var eventCounts = new Dictionary<string, int> {{cnps.Keys.ToList()[0], 2 }};
        var mean = SummaryFeatures.GetMeanBreakpoints(cnps, eventCounts, false);
        Assert.AreEqual(1.5, mean);
    }
    
    [Test]
    public void TestDefaultBreakpointsPerChromosome()
    {
        var karXX = new Karyotype(_genRef, true);
        var cnps = GetCNPs(new List<Karyotype> { karXX });
        var includeSexChromosomes = false;
        var values = SummaryFeatures.GetBreakpointsPerChromosome(_genRef, cnps, includeSexChromosomes);
        // Check that only one sample was found
        Assert.AreEqual(1, values.Count);
        // Check that the number of chromosomes is correct
        Assert.AreEqual(22, values["sample_1"].Count);
        // Check that no breakpoints were found for any of the chromosomes
        foreach (var val in values["sample_1"])
        {
            Assert.AreEqual(0, val);
        }
    }

    [Test]
    public void TestBreakpointsPerChromosome()
    {
        var karXX = new Karyotype(_genRef, true);
        // Missing segment from the start
        karXX.ApplyInternalDeletion(0, 0, 1000);
        var cnps = GetCNPs(new List<Karyotype> { karXX });
        var includeSexChromosomes = false;
        var values = SummaryFeatures.GetBreakpointsPerChromosome(_genRef, cnps, includeSexChromosomes);
        Assert.AreEqual(1, values.Count);
        // Chromosome 1 should have 1 breakpoint from the deletion
        Assert.AreEqual(1, values["sample_1"][0]);
        // All others 0
        var nChrs = values["sample_1"].Count;
        for (int i = 1; i < nChrs; i++)
        {
            Assert.AreEqual(0, values["sample_1"][i]);
        }

        // Whole missing internal segment
        karXX = new Karyotype(_genRef, true);
        // Missing segment from the start
        karXX.ApplyInternalDeletion(0, 1000, 2000);
        cnps = GetCNPs(new List<Karyotype> { karXX });
        values= SummaryFeatures.GetBreakpointsPerChromosome(_genRef, cnps, includeSexChromosomes);
        Assert.AreEqual(1, values.Count);
        // Chromosome 1 should have 2 breakpoints from the deletion
        Assert.AreEqual(2, values["sample_1"][0]);
        // All others 0
        for (int i = 1; i < nChrs; i++)
        {
            Assert.AreEqual(0, values["sample_1"][i]);
        }
    }

    [Test]
    public void TestGetBreakpointDistribution()
    {
        var karXX = new Karyotype(_genRef, true);
        // Create some breakpoints in one sample
        // 2 breakpoints on chr1
        karXX.ApplyInternalDeletion(0, 0, 1000);
        // 4 breakpoints on chr2
        karXX.ApplyInternalDuplication(1, 1000, 2000);
        karXX.ApplyInternalDuplication(1, 5000, 6000);
        var karXY = new Karyotype(_genRef, false);
        var cnps = GetCNPs(new List<Karyotype> { karXX, karXY });
        var includeSexChromosomes = false;
        var values = SummaryFeatures.GetBreakpointsDistribution(_genRef, cnps, includeSexChromosomes);
        // One entry for each autosome
        Assert.AreEqual(22, values.Count);
        Assert.AreEqual(0.5, values[0], double.Epsilon);
        Assert.AreEqual(2, values[1], double.Epsilon);
        for (int i = 2; i < values.Count; i++)
        {
            Assert.AreEqual(0, values[i], double.Epsilon);
        }
    }

    [Test]
    public void TestDefaultBreakpoints()
    {
        var karXX = new Karyotype(_genRef, true);
        var karXY = new Karyotype(_genRef, false);
        var cnps = GetCNPs(new List<Karyotype> { karXX, karXY });
        var size = 10_000_000;
        var includeSexChromosomes = false;
        var values = SummaryFeatures.GetBreakpointsPerBin(_genRef, cnps, includeSexChromosomes, size);
        Assert.AreEqual(2, values.Count);
        // The number of bins across the chromosome 
        Assert.AreEqual(300, values["sample_1"].Count);
        Assert.AreEqual(300, values["sample_2"].Count);
        for (int i = 0; i < 300; i++)
        {
            Assert.AreEqual(0, values["sample_1"][i]);
            Assert.AreEqual(0, values["sample_2"][i]);
        }
        Assert.AreEqual(0, values["sample_1"].Count(x => x == 1));
        Assert.AreEqual(0, values["sample_2"].Count(x => x == 1));
    }

    [Test]
    public void TestBreakpoints()
    {
        var karXX = new Karyotype(_genRef, true);
        // Missing segment from the start
        karXX.ApplyInternalDeletion(0, 0, 1000);
        var cnps = GetCNPs(new List<Karyotype> { karXX});
        var size = 10_000_000;
        var includeSexChromosomes = false;
        var values = SummaryFeatures.GetBreakpointsPerBin(_genRef, cnps, includeSexChromosomes, size);
        Assert.AreEqual(1, values.Count);
        Assert.AreEqual(300, values["sample_1"].Count);
        // First bin should have a breakpoint
        Assert.AreEqual(1, values["sample_1"][0]);
        // All others are empty
        for (int i = 1; i < 300; i++)
        {
            Assert.AreEqual(0, values["sample_1"][i]);
        }
        
        // This time with a totally missing segment
        karXX = new Karyotype(_genRef, true);
        karXX.ApplyInternalDeletion(0, 1000, 2000);
        cnps = GetCNPs(new List<Karyotype> { karXX});
        values = SummaryFeatures.GetBreakpointsPerBin(_genRef, cnps, includeSexChromosomes, size);
        Assert.AreEqual(1, values.Count);
        Assert.AreEqual(300, values["sample_1"].Count);
        // First bin should have two breakpoints (start and end of missing segment)
        Assert.AreEqual(2, values["sample_1"][0]);
        // All others are empty
        for (int i = 1; i < 300; i++)
        {
            Assert.AreEqual(0, values["sample_1"][i]);
        }
    }

    [Test]
    public void TestDefaultMajMinCNs()
    {
        var karXX = new Karyotype(_genRef, true);
        var cnps = GetCNPs(new List<Karyotype> { karXX });
        // Major copy number
        var (values, max) = SummaryFeatures.GetMajMinCNs(cnps, true);
        // One sample
        Assert.AreEqual(1, values.Count);
        // Average is 1
        Assert.AreEqual(1, values[0]);
        // Maximum is 1
        Assert.AreEqual(1, max);
        // Minor copy number
        (values, max) = SummaryFeatures.GetMajMinCNs(cnps, false);
        Assert.AreEqual(1, values.Count);
        Assert.AreEqual(1, values[0]);
        Assert.AreEqual(1, max);
    }

    [Test]
    public void TestMajMinCNs()
    {
        // Create a karyotype with 2 copies of chr1-h1, 
        // and no other chromosomes
        var karXX = new Karyotype(_genRef, true);
        karXX.ApplyContigDuplication(0);
        for (int i = 1; i < 46; i++)
        {
            if (i == 23)
            {
                continue;
            }
            karXX.ApplyContigDeletion(i);
        }
        var cnps = GetCNPs(new List<Karyotype> { karXX });
        // Major copy number
        var (values, max) = SummaryFeatures.GetMajMinCNs(cnps, true);
        Assert.AreEqual(1, values.Count);
        // 2x chr1-h1, 1x chr2-h1
        var chr1Length = _genRef.ChrLengths[_genRef.AllChrs[0]];
        var haploidLength = _genRef.GetGenomeLen(karXX.SexXX,false);
        var majFrac = 2.0*chr1Length/haploidLength;
        Assert.AreEqual(majFrac, values[0], double.Epsilon);
        Assert.AreEqual(majFrac, max, double.Epsilon);
        // Minor copy number
        (values, max) = SummaryFeatures.GetMajMinCNs(cnps, false);
        var minFrac = chr1Length/(double)haploidLength;
        Assert.AreEqual(1, values.Count);
        Assert.AreEqual(minFrac, values[0], double.Epsilon);
        Assert.AreEqual(minFrac, max, double.Epsilon);
    }

    [Test]
    public void TestDefaultHomozygousDeletionFraction()
    {
        var karXX = new Karyotype(_genRef, true);
        var cnps = GetCNPs(new List<Karyotype> { karXX });
        var values = SummaryFeatures.GetHomozygousDeletionFraction(_genRef, cnps);
        // One sample
        Assert.AreEqual(1, values.Count);
        // No homozygous deletions
        Assert.AreEqual(0, values[0], double.Epsilon);
    }
    [Test]
    public void TestHomozygousDeletionFraction()
    {
        var karXX = new Karyotype(_genRef, true);
        karXX.ApplyContigDeletion(0);
        karXX.ApplyContigDeletion(23);
        var chr1Length = _genRef.ChrLengths[_genRef.AllChrs[0]];
        var autosomeLen = _genRef.AutosomeLen;
        var homoDelFrac = chr1Length/(double)autosomeLen;
        var cnps = GetCNPs(new List<Karyotype> { karXX });
        var values = SummaryFeatures.GetHomozygousDeletionFraction(_genRef, cnps);
        Assert.AreEqual(1, values.Count);
        Assert.AreEqual(homoDelFrac, values[0], double.Epsilon);
    }
    [Test]
    public void TestDefaultMeanCNAlongGenome()
    {
        var karXX = new Karyotype(_genRef, true);
        // Bin the genome into 1Mb bins
        var binWidth = 1_000_000;
        var includeSexChromosomes = false;
        var binner = new Binner(_genRef, binWidth, includeSexChromosomes);
        var karDict = new Dictionary<string, Karyotype> {{"sample_1", karXX}};
        var binnedCNPs = binner.GetBinnedCNProfiles(karDict);
        //IsFemaleSimulatedDict = samples.ToDictionary(s => s.SampleId, s => s.SexXX);
        var values = SummaryFeatures.GetMeanCNAlongGenome(binnedCNPs);
        // Without sex chromosomes
        Assert.AreEqual(2897, values.Count);
        foreach (var val in values)
        {
            Assert.AreEqual(2, val);
        }
        // with sex chromosomes
        includeSexChromosomes = true;
        binner = new Binner(_genRef, binWidth, includeSexChromosomes);
        binnedCNPs = binner.GetBinnedCNProfiles(karDict);
        values = SummaryFeatures.GetMeanCNAlongGenome(binnedCNPs);
        Assert.AreEqual(3113, values.Count);
        for (int i = 0; i < 2897; i++)
        {
            Assert.AreEqual(2, values[i]);
        }
        // X chromosome
        for (int i = 0; i < 156; i++)
        {
            Assert.AreEqual(2, values[2897+i]);
        }
        // Y chromosome
        for (int i = 0; i < 60; i++)
        {
            Assert.AreEqual(0, values[3053+i]);
        }
    }

    [Test]
    public void TestMeanCNAlongGenome()
    {
        var karXX = new Karyotype(_genRef, true);
        karXX.ApplyContigDuplication(0);
        karXX.ApplyContigDeletion(1);
        karXX.ApplyInternalDeletion(2, 1000, 200000);
        // Bin the genome into 1Mb bins
        var binWidth = 1_000_000;
        var includeSexChromosomes = false;
        var binner = new Binner(_genRef, binWidth, includeSexChromosomes);
        var karDict = new Dictionary<string, Karyotype> {{"sample_1", karXX}};
        var binnedCNPs = binner.GetBinnedCNProfiles(karDict);
        var values = SummaryFeatures.GetMeanCNAlongGenome(binnedCNPs);
        Assert.AreEqual(2897, values.Count);
        // Chrom Duplication on chr1
        for (int i = 0; i < 250; i++)
        {
            Assert.AreEqual(3, values[i]);
        }
        // Chrom deletion on chr2
        for (int i = 250; i < 250+244; i++)
        {
            Assert.AreEqual(1, values[i]);
        }
        var n = values.GetRange(494, 494+199);
        // Segmental deletion on chr3
        for (int i = 494; i < 494+2; i++)
        {
            Assert.AreEqual(1, values[i]);
        }
        // Rest of chr3 is cn-normal
        for (int i = 496; i < 494+199; i++)
        {
            Assert.AreEqual(2, values[i]);
        }
        // Rest of the autosomes
        for (int i = 693; i < 2897; i++)
        {
            Assert.AreEqual(2, values[i]);
        }
    }
}