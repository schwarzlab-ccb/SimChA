using System;
using System.Collections.Generic;
using NUnit.Framework;
using SimChA.Computation;
using SimChA.Data;

namespace Tests;

[TestFixture]
public class TestRegions
{
    private Region _cRegion;

    [SetUp]
    public void Setup()
    {
        _cRegion = new Region(0, 249250621, "chr1", true, [], []);
    }

    [Test]
    public void TestOverlap()
    {
        var testRange = new GenRange(1000, 2000, "chr1");
        Assert.IsTrue(testRange.Overlaps(new GenRange(500, 1500, "chr1")));
        Assert.IsFalse(testRange.Overlaps(new GenRange(500, 1500, "chr2")));
        
        Assert.IsTrue(testRange.Overlaps(new GenRange(0, 1001, "chr1")));
        Assert.IsFalse(testRange.Overlaps(new GenRange(0, 1000, "chr1")));

        Assert.IsTrue(testRange.Overlaps(new GenRange(1999, 3000, "chr1")));
        Assert.IsFalse(testRange.Overlaps(new GenRange(2000, 3000, "chr1")));

        var gene = new Gene(69090, 70008, "chr1" , GeneLT.TSG, 0, 0.142321064 );
        var range = new GenRange(0, 249250621, "chr1");
        Assert.IsTrue(range.Overlaps(gene));
    }

    [Test]
    public void TestDeletion()
    {
        var regions = new List<Region> { _cRegion };
        // Start and end within a region
        regions = RegionOps.DeleteRange(regions, 1000, 2000);
        Console.WriteLine(Contig.ToString(regions));
        Assert.AreEqual(_cRegion.Length - 1000, RegionOps.CountLength(regions));
        // Start and end within a region (out of two)
        regions = RegionOps.DeleteRange(regions, 2000, 4000);
        Console.WriteLine(Contig.ToString(regions));
        Assert.AreEqual(_cRegion.Length - 3000, RegionOps.CountLength(regions));
        // Two neighbouring regions
        regions = RegionOps.DeleteRange(regions, 1500, 2500);
        Console.WriteLine(Contig.ToString(regions));
        Assert.AreEqual(_cRegion.Length - 4000, RegionOps.CountLength(regions));
        // Cut region out
        regions = RegionOps.DeleteRange(regions, 500, 2500);
        Console.WriteLine(Contig.ToString(regions));
        Assert.AreEqual(_cRegion.Length - 6000, RegionOps.CountLength(regions));
        // Remove region from front
        regions = RegionOps.DeleteRange(regions, 0, 500);
        Console.WriteLine(Contig.ToString(regions));
        Assert.AreEqual(_cRegion.Length - 6500, RegionOps.CountLength(regions));
        // oversized range selection
        regions = RegionOps.DeleteRange(regions, -1000, _cRegion.Length);
        Console.WriteLine(Contig.ToString(regions));
        Assert.AreEqual(0, RegionOps.CountLength(regions));
    }
    
    [Test]
    public void TestDeleteInverted()
    {
        var regions = new List<Region> { _cRegion };
        RegionOps.Revert(regions);
        var newRegions = RegionOps.DeleteRange(regions, 1000, 2000);
        Assert.AreEqual(1000, newRegions[0].Length);
        Assert.AreEqual(_cRegion.Length - 2000, newRegions[1].Length);
        Assert.AreEqual(-_cRegion.Length, newRegions[0].Start);
        Assert.AreEqual(-_cRegion.Length + 1000, newRegions[0].End);
        Assert.AreEqual(-_cRegion.Length + 2000, newRegions[1].Start);
        Assert.AreEqual(0, newRegions[1].End);
    }

    [Test]
    public void TestDeletionUnitLength()
    {
        var regions = new List<Region>
        {

            new(0, 1, "chr1", true, [], []),
            new(1, 2, "chr1", true, [], []),
            new(2, 3, "chr1", true, [], []),
            new(3, 4, "chr1", true, [], [])
        };
        var res = new List<Region>
        {
            new(0, 1, "chr1", true, [], []),
            new(3, 4, "chr1", true, [], [])
        };
        Assert.AreEqual(res, RegionOps.DeleteRange(regions, 1, 3));
    }
    
    [Test]
    public void TestMergeRegions()
    {
        List<Region> regions =
        [
            new(0, 1, "chr1", true, [], []),
            new(1, 2, "chr1", true, [], []),
            new(2, 3, "chr1", true, [], []),
            new(3, 4, "chr1", true, [], [])
        ];
        List<Region> mergeAll =
        [
            new(0, 4, "chr1", true, [], [])
        ];
        Assert.AreEqual(mergeAll, RegionOps.MergeRegions(regions));
        
        List<Region> mergeAllInverted =
        [
            new(-4, -2, "chr1", true, [], []),
            new(-2, 0, "chr1", true, [], [])
        ];
        var merged = RegionOps.MergeRegions(mergeAllInverted);
        RegionOps.Revert(merged);
        Assert.AreEqual(mergeAll, merged);
        mergeAllInverted =
        [
            new(-4, -2, "chr1", true, [], []),
            new(-2, 0, "chr1", true, [], [])
        ];
        RegionOps.Revert(mergeAllInverted);
        merged =  RegionOps.MergeRegions(mergeAllInverted);
        Assert.AreEqual(mergeAll, merged);

        regions =
        [
            new(0, 1, "chr1", true, [], []),
            new(1, 2, "chr1", true, [], []),
            new(-3, -2, "chr1", true, [], []),
            new(3, 4, "chr1", true, [], [])
        ];
        List<Region> mergeOneInverted =
        [
            new(0, 2, "chr1", true, [], []),
            new(-3, -2, "chr1", true, [], []),
            new(3, 4, "chr1", true, [], [])
        ];
        Assert.AreEqual(mergeOneInverted, RegionOps.MergeRegions(regions));

        regions =
        [
            new(0, 1, "chr1", true, [], []),
            new(1, 2, "chr1", true, [], []),
            new(3, 4, "chr1", true, [], [])
        ];
        List<Region> mergeGapped =
        [
            new(0, 2, "chr1", true, [], []),
            new(3, 4, "chr1", true, [], [])
        ];
        Assert.AreEqual(mergeGapped, RegionOps.MergeRegions(regions));
    }

    [Test]
    public void TestCopy()
    {
        var regions = new List<Region> { _cRegion };
        // Start and end within a region
        var regCopy = RegionOps.CopyRange(regions, 500, 3500);
        Console.WriteLine(Contig.ToString(regCopy));
        Assert.AreEqual(3000, RegionOps.CountLength(regCopy));
        // Copy across regions
        regions = RegionOps.DeleteRange(regions, 1000, 2000);
        regions = RegionOps.DeleteRange(regions, 4000, 5000);
        regCopy = RegionOps.CopyRange(regions, 500, 3500);
        Console.WriteLine(Contig.ToString(regCopy));
        Assert.AreEqual(3000, RegionOps.CountLength(regCopy));
    }

    [Test]
    public void TestCopyUnitLength()
    {
        var regions = new List<Region>
        {
            new(0, 1, "chr1", true, [], []),
            new(1, 2, "chr1", true, [], []),
            new(2, 3, "chr1", true, [], []),
            new(3, 4, "chr1", true, [], [])
        };
        var res = new List<Region>
        {
            new(1, 2, "chr1", true, [], []),
            new(2, 3, "chr1", true, [], []),
        };
        var copy = RegionOps.CopyRange(regions, 1, 3);
        Assert.AreEqual(res, copy);
    }

    [Test]
    public void TestCopyInverted()
    {
        var regions = new List<Region> { _cRegion };
        RegionOps.Revert(regions);
        // Start and end within a region
        var regCopy = RegionOps.CopyRange(regions, 1000, 2000);
        Assert.AreEqual(1, regions.Count);
        Assert.AreEqual(1000, regCopy[0].Length);
        Assert.AreEqual(-_cRegion.Length + 1000, regCopy[0].Start);
        Assert.AreEqual(-_cRegion.Length + 2000, regCopy[0].End);
    }

    [Test]
    public void TestSplit()
    {
        var regions = new List<Region> { _cRegion };
        var (before, after) = RegionOps.SplitRegions(regions, 2000);
        Assert.AreEqual(2000, before[0].Length);
        Assert.AreEqual(2000, before[0].End);
        Assert.AreEqual(2000, after[0].Start);

        RegionOps.Revert(regions);
        (before, after) = RegionOps.SplitRegions(regions, 2000);
        Assert.AreEqual(2000, before[0].Length);
        Assert.AreEqual(_cRegion.Length - 2000,  after[0].Length);
    }
    
    [Test]
    public void TestSplitUnitLength()
    {
        var regions = new List<Region>
        {
            new(0, 2, "chr1", true, [], [])
        };
        var (before, after) = RegionOps.SplitRegions(regions, 1);
        Assert.AreEqual(1, RegionOps.CountLength(after));
    }
    
    [Test]
    public void TestInverse()
    {
        var regions = new List<Region> { _cRegion };
        regions = RegionOps.DeleteRange(regions, 1000, 2000);
        RegionOps.Revert(regions);
        Console.WriteLine(Contig.ToString(regions));
        Assert.AreEqual(_cRegion.Length - 1000, RegionOps.CountLength(regions));
    }
    
    [Test]
    public void TestConcat()
    {
        var regions = new List<Region> { _cRegion };
        regions = RegionOps.ConcatRegions(regions, regions);
        Console.WriteLine(Contig.ToString(regions));
        Assert.AreEqual(_cRegion.Length * 2, RegionOps.CountLength(regions));
        regions = RegionOps.ConcatRegions(new[] {regions, regions });
        Console.WriteLine(Contig.ToString(regions));
        Assert.AreEqual(_cRegion.Length * 4, RegionOps.CountLength(regions));
    }

    [Test]
    public void TestPointMutateRegion()
    {
        var regions = new List<Region>
        {
            new(0, 1, "chr1", true, [], []),
            new(1, 2, "chr1", true, [], []),
            new(2, 3, "chr1", true, [], []),
            new(3, 4, "chr1", true, [], [])
        };
        const int location = 2;
        const Nucleotide newNucleotide = Nucleotide.A;
        const Nucleotide oldNucleotide = Nucleotide.C;
        var mutatedRegions = RegionOps.Copy(regions);
        RegionOps.PointMutateRegion(mutatedRegions, location, oldNucleotide, newNucleotide);
        Assert.AreEqual(regions[0], mutatedRegions[0]);
        Assert.AreEqual(regions[1], mutatedRegions[1]);
        Assert.AreNotEqual(regions[2], mutatedRegions[2]);
        Assert.AreEqual(regions[3], mutatedRegions[3]);

        Assert.IsNotNull(mutatedRegions[2].SNVs);
        Assert.AreEqual(1, mutatedRegions[2].SNVs!.Count);
        Assert.AreEqual(location, mutatedRegions[2].SNVs![0].Pos);
        Assert.AreEqual(newNucleotide, mutatedRegions[2].SNVs![0].Alt);

    }

    [Test]
    public void TestUpdateSNVs()
    {
        var regions = new List<Region>
        {
            new(0, 100, "chr1", true, [], []),
        };
        const int location = 45;
        const Nucleotide newNucleotide = Nucleotide.A;
        const Nucleotide oldNucleotide = Nucleotide.C;
        var mutatedRegions = RegionOps.Copy(regions);
        RegionOps.PointMutateRegion(mutatedRegions, location, oldNucleotide, newNucleotide);
        foreach (var region in mutatedRegions)
        {
            Assert.IsNotNull(region.SNVs);
        }
        var finalRegions = RegionOps.DeleteRange(mutatedRegions, 30, 70);
        foreach (var region in finalRegions)
        {
            Assert.IsEmpty(region.SNVs);
        }
    }

    [Test]
    public void TestEquality()
    {
        var region1 = new Region(0, 100, "chr1", true, [], []);
        var region2 = new Region(0, 100, "chr1", true, [], []);
        Assert.AreEqual(region1, region2);
        
        var range1 = new GenRange(0, 100, "chr1");
        var range2 = new GenRange(0, 100, "chr1");
        Assert.AreEqual(range1, range2);
    }
    
    [Test]
    public void TestPresentGenesWithDeletion()
    {
        var tsg = new Gene(0, 10, "chr1", GeneLT.TSG, 0, 0.1);
        var og = new Gene(20, 30, "chr1", GeneLT.OG, 0, 0.2);
        var ess = new Gene(40, 50, "chr1", GeneLT.Ess, 0, 0.3);
        var r = new Region(0, 249250621, "chr1", true, [], [tsg, og, ess]);
        var regions = new List<Region>{r};
        // Delete a region affecting the og
        regions = RegionOps.DeleteRange(regions, 15, 25);
        var tsgRegion = regions[0];
        var essRegion = regions[1];
        // Check that the original region is unchanged
        Assert.AreEqual(1, r.CountGeneType(GeneLT.TSG));
        Assert.AreEqual(1, r.CountGeneType(GeneLT.OG));
        Assert.AreEqual(1, r.CountGeneType(GeneLT.Ess));
        
        // Check that the tsg and essential genes are intact
        Assert.AreEqual(1, tsgRegion.CountGeneType(GeneLT.TSG));
        Assert.AreEqual(0, tsgRegion.CountGeneType(GeneLT.OG));
        Assert.AreEqual(0, tsgRegion.CountGeneType(GeneLT.Ess));

        Assert.AreEqual(0, essRegion.CountGeneType(GeneLT.TSG));
        Assert.AreEqual(0, essRegion.CountGeneType(GeneLT.OG));
        Assert.AreEqual(1, essRegion.CountGeneType(GeneLT.Ess));
    }
}
