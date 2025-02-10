using System;
using System.Collections.Generic;
using Extreme.DataAnalysis.Models;
using NUnit.Framework;
using SimChA.Computation;
using SimChA.Data;
using SimChA.DataTypes;
using SimChA.Simulation;

namespace Tests;

[TestFixture]
public class TestRegions
{
    private Region _cRegion;

    [SetUp]
    public void Setup()
    {
        _cRegion = new Region(0, 249250621, "chr1", true, true);
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

        var gene = new Gene("OR4F5", new GenRange(69090, 70008, "chr1"), 0.142321064);
        var range = new GenRange(0, 249250621, "chr1");
        Assert.IsTrue(range.Overlaps(gene.Range));
    }

    [Test]
    public void TestDeletion()
    {
        var regions = new List<Region> { _cRegion };
        // Start and end within a region
        regions = RegionOps.DeleteRange(regions, 1000, 2000);
        Console.WriteLine(Contig.ToString(regions));
        Assert.AreEqual(_cRegion.Length - 1000, Contig.Length(regions));
        // Start and end within a region (out of two)
        regions = RegionOps.DeleteRange(regions, 2000, 4000);
        Console.WriteLine(Contig.ToString(regions));
        Assert.AreEqual(_cRegion.Length - 3000, Contig.Length(regions));
        // Two neighbouring regions
        regions = RegionOps.DeleteRange(regions, 1500, 2500);
        Console.WriteLine(Contig.ToString(regions));
        Assert.AreEqual(_cRegion.Length - 4000, Contig.Length(regions));
        // Cut region out
        regions = RegionOps.DeleteRange(regions, 500, 2500);
        Console.WriteLine(Contig.ToString(regions));
        Assert.AreEqual(_cRegion.Length - 6000, Contig.Length(regions));
        // Remove region from front
        regions = RegionOps.DeleteRange(regions, 0, 500);
        Console.WriteLine(Contig.ToString(regions));
        Assert.AreEqual(_cRegion.Length - 6500, Contig.Length(regions));
        // oversized range selection
        regions = RegionOps.DeleteRange(regions, -1000, _cRegion.Length);
        Console.WriteLine(Contig.ToString(regions));
        Assert.AreEqual(0, Contig.Length(regions));
    }
    
    [Test]
    public void TestDeleteInverted()
    {
        var regions = RegionOps.InvertRegions(new List<Region> { _cRegion });
        var newRegions = RegionOps.DeleteRange(regions, 1000, 2000);
        Assert.AreEqual(1000, newRegions[0].Length);
        Assert.AreEqual(_cRegion.Length - 2000, newRegions[1].Length);
        Assert.AreEqual(_cRegion.Length - 1000, newRegions[0].Start);
        Assert.AreEqual(_cRegion.Length, newRegions[0].End);
        Assert.AreEqual(0, newRegions[1].Start);
        Assert.AreEqual(_cRegion.Length - 2000, newRegions[1].End);
    }

    [Test]
    public void TestDeletionUnitLength()
    {
        var regions = new List<Region>
        {
            _cRegion with { Start = 0, End = 1 },
            _cRegion with { Start = 1, End = 2 },
            _cRegion with { Start = 2, End = 3 },
            _cRegion with { Start = 3, End = 4 }
        };
        var res = new List<Region>
        {
            _cRegion with { Start = 0, End = 1 },
            _cRegion with { Start = 3, End = 4 }
        };
        Assert.AreEqual(res, RegionOps.DeleteRange(regions, 1, 3));
    }
    
    [Test]
    public void TestDeleteArm()
    {
        var regions = new List<Region>
        {
            new PArm(0, 1, "chr1", true, true),
            new Centromere(1, 2, "chr1", true, true),
            new QArm(2, 3, "chr1", true, true)
        };
        var pArm = new List<Region>
        {
            new Centromere(1, 2, "chr1", true, true),
            new QArm(2, 3, "chr1", true, true)
        };
        Assert.AreEqual(pArm, RegionOps.DeleteArm(regions, 1, true, false));
        regions = new List<Region>
        {
            new PArm(0, 1, "chr1", true, true),
            new Centromere(1, 2, "chr1", true, true),
            new QArm(2, 3, "chr1", true, true)
        };
        var qArm = new List<Region>
        {
            new PArm(0, 1, "chr1", true, true),
            new Centromere(1, 2, "chr1", true, true)
        };
        Assert.AreEqual(qArm, RegionOps.DeleteArm(regions, 1, false, false));
        regions = new List<Region>
        {
            new PArm(0, 1, "chr1", true, true),
            new Centromere(1, 2, "chr1", true, true),
            new QArm(2, 3, "chr1", true, true)
        };
        var pArmWithCentromere = new List<Region>
        {
            new QArm(2, 3, "chr1", true, true)
        };
        Assert.AreEqual(pArmWithCentromere, RegionOps.DeleteArm(regions, 1, true, true));
        regions = new List<Region>
        {
            new PArm(0, 1, "chr1", true, true),
            new Centromere(1, 2, "chr1", true, true),
            new QArm(2, 3, "chr1", true, true)
        };
        var qArmWithCentromere = new List<Region>
        {
            new PArm(0, 1, "chr1", true, true)
        };
        Assert.AreEqual(qArmWithCentromere, RegionOps.DeleteArm(regions, 1, false, true));
    }

    [Test]
    public void TestGetArm()
    {
        var regions = new List<Region>
        {
            new PArm(0, 1, "chr1", true, true),
            new Centromere(1, 2, "chr1", true, true),
            new QArm(2, 3, "chr1", true, true)
        };
        var qArm = new List<Region>
        {
            new QArm(2, 3, "chr1", true, true)
        };
        Assert.AreEqual(qArm, RegionOps.GetArm(regions, 1, false, false));
        var pArm = new List<Region>
        {
            new PArm(0, 1, "chr1", true, true)
        };
        Assert.AreEqual(pArm, RegionOps.GetArm(regions, 1, true, false));
        var qArmWithCentromere = new List<Region>
        {
            new Centromere(1, 2, "chr1", true, true),
            new QArm(2, 3, "chr1", true, true)
        };
        Assert.AreEqual(qArmWithCentromere, RegionOps.GetArm(regions, 1, false, true));
        var pArmWithCentromere = new List<Region>
        {
            new PArm(0, 1, "chr1", true, true),
            new Centromere(1, 2, "chr1", true, true)
        };
        Assert.AreEqual(pArmWithCentromere, RegionOps.GetArm(regions, 1, true, true));        
    }

    [Test]
    public void TestMergeRegions()
    {
        var regions = new List<Region>
        {
            _cRegion with { Start = 0, End = 1 },
            _cRegion with { Start = 1, End = 2 },
            _cRegion with { Start = 2, End = 3 },
            _cRegion with { Start = 3, End = 4 }
        };
        var mergeAll = new List<Region>
        {
            _cRegion with { Start = 0, End = 4 },
        };
        Assert.AreEqual(mergeAll, RegionOps.MergeRegions(regions));

        regions = new List<Region>
        {
            _cRegion with { Start = 0, End = 1 },
            _cRegion with { Start = 1, End = 2 },
            _cRegion with { Start = 2, End = 3, Forward = false },
            _cRegion with { Start = 3, End = 4, Forward = false }
        };
        var mergeAllInverted = new List<Region>
        {
            _cRegion with { Start = 0, End = 2},
            _cRegion with { Start = 2, End = 4, Forward = false },
        };
        Assert.AreEqual(mergeAllInverted, RegionOps.MergeRegions(regions));

        regions = new List<Region>
        {
            _cRegion with { Start = 0, End = 1 },
            _cRegion with { Start = 1, End = 2 },
            _cRegion with { Start = 2, End = 3, Forward = false },
            _cRegion with { Start = 3, End = 4}
        };
        var mergeOneInverted = new List<Region>
        {
            _cRegion with { Start = 0, End = 2},
            _cRegion with { Start = 2, End = 3, Forward = false },
            _cRegion with { Start = 3, End = 4}
        };
        Assert.AreEqual(mergeOneInverted, RegionOps.MergeRegions(regions));

        regions = new List<Region>
        {
            _cRegion with { Start = 0, End = 1 },
            _cRegion with { Start = 1, End = 2 },
            _cRegion with { Start = 3, End = 4}
        };
        var mergeGapped = new List<Region>
        {
            _cRegion with { Start = 0, End = 2},
            _cRegion with { Start = 3, End = 4}
        };
        Assert.AreEqual(mergeGapped, RegionOps.MergeRegions(regions));
    }

    [Test]
    public void TestCopy()
    {
        var regions = new List<Region> { _cRegion };
        // Start and end within a region
        var regCopy = RegionOps.CopyRange(regions, 500, 3500);
        Console.WriteLine(Contig.ToString(regCopy));
        Assert.AreEqual(3000, Contig.Length(regCopy));
        // Copy across regions
        regions = RegionOps.DeleteRange(regions, 1000, 2000);
        regions = RegionOps.DeleteRange(regions, 4000, 5000);
        regCopy = RegionOps.CopyRange(regions, 500, 3500);
        Console.WriteLine(Contig.ToString(regCopy));
        Assert.AreEqual(3000, Contig.Length(regCopy));
    }

    [Test]
    public void TestCopyUnitLength()
    {
        var regions = new List<Region>
        {
            _cRegion with { Start = 0, End = 1 },
            _cRegion with { Start = 1, End = 2 },
            _cRegion with { Start = 2, End = 3 },
            _cRegion with { Start = 3, End = 4 }
        };
        var res = new List<Region>
        {
            _cRegion with { Start = 1, End = 2 },
            _cRegion with { Start = 2, End = 3 }
        };
        Assert.AreEqual(res, RegionOps.CopyRange(regions, 1, 3));
    }

    [Test]
    public void TestCopyInverted()
    {
        var regions = RegionOps.InvertRegions(new List<Region> { _cRegion });
        // Start and end within a region
        var regCopy = RegionOps.CopyRange(regions, 1000, 2000);
        Assert.AreEqual(1, regions.Count);
        Assert.AreEqual(1000, regCopy[0].Length);
        Assert.AreEqual(_cRegion.Length - 2000, regCopy[0].Start);
        Assert.AreEqual(_cRegion.Length - 1000, regCopy[0].End);
    }

    [Test]
    public void TestSplit()
    {
        var regions = new List<Region> { _cRegion };
        
        var (before, after) = RegionOps.SplitRegions(regions, 2000);
        Assert.AreEqual(2000, before[0].Length);
        Assert.AreEqual(regions[0]  with { End = 2000}, before[0]);
        Assert.AreEqual(regions[0]  with { Start = 2000}, after[0]);

        regions = RegionOps.InvertRegions(regions);
        (before, after) = RegionOps.SplitRegions(regions, 2000);
        Assert.AreEqual(2000, before[0].Length);
        Assert.AreEqual(regions[0] with { Start = _cRegion.Length - 2000}, before[0]);
        Assert.AreEqual(regions[0] with { End = _cRegion.Length - 2000}, after[0]);
    }
    
    [Test]
    public void TestSplitUnitLength()
    {
        var regions = new List<Region>
        {
            _cRegion with { Start = 0, End = 2 }
        };
        var (before, after) = RegionOps.SplitRegions(regions, 1);
        Assert.AreEqual(1, Contig.Length(after));
    }
    
    [Test]
    public void TestInverse()
    {
        var regions = new List<Region> { _cRegion };
        regions = RegionOps.DeleteRange(regions, 1000, 2000);
        regions = RegionOps.InvertRegions(regions);
        Console.WriteLine(Contig.ToString(regions));
        Assert.AreEqual(_cRegion.Length - 1000, Contig.Length(regions));
    }
    
    [Test]
    public void TestConcat()
    {
        var regions = new List<Region> { _cRegion };
        regions = RegionOps.ConcatRegions(regions, regions);
        Console.WriteLine(Contig.ToString(regions));
        Assert.AreEqual(_cRegion.Length * 2, Contig.Length(regions));
        regions = RegionOps.ConcatRegions(new[] {regions, regions });
        Console.WriteLine(Contig.ToString(regions));
        Assert.AreEqual(_cRegion.Length * 4, Contig.Length(regions));
    }

    [Test]
    public void TestPointMutateRegion()
    {
        var regions = new List<Region>
        {
            _cRegion with { Start = 0, End = 1 },
            _cRegion with { Start = 1, End = 2 },
            _cRegion with { Start = 2, End = 3 },
            _cRegion with { Start = 3, End = 4 }
        };
        var location = 2;
        var newNucleotide = Nucleotide.A;
        var mutatedRegions = RegionOps.PointMutateRegion(regions, location, newNucleotide);
        Assert.AreEqual(regions[0], mutatedRegions[0]);
        Assert.AreEqual(regions[1], mutatedRegions[1]);
        Assert.AreNotEqual(regions[2], mutatedRegions[2]);
        Assert.AreEqual(regions[3], mutatedRegions[3]);

        Assert.NotNull(mutatedRegions[2].SNVDict);
        Assert.AreEqual(newNucleotide, mutatedRegions[2].SNVDict[location]);

    }

    [Test]
    public void TestFindRegion()
    {
        var regions = new List<Region>
        {
            _cRegion with { Start = 0, End = 1 },
            _cRegion with { Start = 1, End = 2 },
            _cRegion with { Start = 2, End = 3 },
            _cRegion with { Start = 3, End = 4 }
        };
        var location = 2;
        (var region, var internalLocation) = RegionOps.FindRegion(regions, location);

        Assert.AreEqual(regions[2], region);
        Assert.AreEqual(location, internalLocation);

    }

    [Test]
    public void TestUpdateSNVDict()
    {
        var regions = new List<Region>
        {
            _cRegion with {Start = 0, End = 100}
        };
        var location = 45;
        var newNucleotide = Nucleotide.A;
        var mutatedRegions = RegionOps.PointMutateRegion(regions, location, newNucleotide);
        foreach (var region in mutatedRegions)
        {
            Assert.IsNotNull(region.SNVDict);
        }
        // TODO: Why does DeleteRange remove the original mutatedRegions SNVDict?
        // Does it matter?
        var finalRegions = RegionOps.DeleteRange(mutatedRegions, 30, 70);
        foreach (var region in finalRegions)
        {
            Assert.IsNull(region.SNVDict);
        }
    }
}
