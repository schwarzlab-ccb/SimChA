// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SimChA.Computation;
using SimChA.Data;
using SimChA.IO;

namespace Tests;

[TestFixture]
public class TestFitness
{
    private const double EPSILON = 0.00001;

    private List<RefGen> _refs = null!;

    [SetUp]
    public void Setup()
    {
        _refs = [
            FileIO.ReadGenRef(TestParsing.DATA_PATH, TestParsing.HG_19, TestParsing.GENE_SET), 
            FileIO.ReadGenRef(TestParsing.DATA_PATH, TestParsing.HG_38, TestParsing.GENE_SET)
        ];
    }
    
    private static Gene MakeGene(string chrNo, double deltaFitness, int index)
        => new(0, 50, chrNo, GeneLT.OG, index, deltaFitness);

    [Test]
    public void TestEssTerm([Values] SexType sex, [Values(0,1)] int refId)
    {
        // Reference copy number per gene, parallel to the gene list: 2 on an autosome.
        int[] diploid = [2];
        Assert.AreEqual(0, Fitness.EssTerm([], [], []), EPSILON);

        // A zero-scored gene contributes nothing whatever its copy number.
        Gene[] testNoEffect = [MakeGene("chr1", 0, 0)];
        Assert.AreEqual(0, Fitness.EssTerm(testNoEffect, [0], diploid), EPSILON);

        // Below the reference dosage costs the gene's score; at or above it costs nothing, and
        // gains are not rewarded.
        Gene[] testGene = [MakeGene("chr1", 0.1, 0)];
        Assert.AreEqual(-0.1, Fitness.EssTerm(testGene, [1], diploid), EPSILON);
        Assert.AreEqual(0, Fitness.EssTerm(testGene, [2], diploid), EPSILON);
        Assert.AreEqual(0, Fitness.EssTerm(testGene, [5], diploid), EPSILON);

        // CN 0 is priced the same as CN 1, not skipped: EvoSimulator forbids it, but the modes that
        // do not must never find losing both copies cheaper than losing one.
        Assert.AreEqual(-0.1, Fitness.EssTerm(testGene, [0], diploid), EPSILON);

        // A hemizygous reference -- a male's X -- is not haploinsufficient at one copy. Scoring
        // against a flat 2 charged every male for his own baseline, which is what the Male cases of
        // TestCalculate caught.
        int[] hemizygous = [1];
        Assert.AreEqual(0, Fitness.EssTerm(testGene, [1], hemizygous), EPSILON);
        Assert.AreEqual(-0.1, Fitness.EssTerm(testGene, [0], hemizygous), EPSILON);

        Gene[] testList = [
            MakeGene("chr1", 0.1, 0),
            MakeGene("chr2", 0.2, 1)
        ];
        int[] bothDiploid = [2, 2];
        Assert.AreEqual(-0.1 + -0.2, Fitness.EssTerm(testList, [1, 1], bothDiploid), EPSILON);
        Assert.AreEqual(-0.2, Fitness.EssTerm(testList, [2, 1], bothDiploid), EPSILON);
        Assert.AreEqual(0, Fitness.EssTerm(testList, [2, 2], bothDiploid), EPSILON);
        // Mixed references: the autosomal gene is short, the hemizygous one is at its baseline.
        Assert.AreEqual(-0.1, Fitness.EssTerm(testList, [1, 1], [2, 1]), EPSILON);
    }

    [Test]
    public void TestAnyEssentialLost()
    {
        Assert.IsFalse(Fitness.AnyEssentialLost([], []));

        Gene[] one = [MakeGene("chr1", 0.1, 0)];
        Assert.IsTrue(Fitness.AnyEssentialLost(one, [0]));
        Assert.IsFalse(Fitness.AnyEssentialLost(one, [1]));
        Assert.IsFalse(Fitness.AnyEssentialLost(one, [2]));

        // Detects a single loss anywhere in the list, and a zero score does not exempt a gene --
        // the prohibition is about viability, not about how much the gene is worth.
        Gene[] many = [
            MakeGene("chr1", 0.1, 0),
            MakeGene("chr2", 0.0, 1),
            MakeGene("chr3", 0.2, 2)
        ];
        Assert.IsFalse(Fitness.AnyEssentialLost(many, [2, 2, 2]));
        Assert.IsTrue(Fitness.AnyEssentialLost(many, [2, 2, 0]));
        Assert.IsTrue(Fitness.AnyEssentialLost(many, [2, 0, 2]));
        Assert.IsTrue(Fitness.AnyEssentialLost(many, [0, 1, 2]));
    }

    [Test]
    public void TestZygosity()
    {
        // Hemizygous
        Assert.AreEqual(0, Fitness.Zygosity([], [], 1), EPSILON);
        // Nullizygous
        Assert.AreEqual(0, Fitness.Zygosity([], [], 0), EPSILON);
        
        Gene[] testNoEffect = [MakeGene("chr1", 0, 0)];
        Assert.AreEqual(0, Fitness.Zygosity(testNoEffect, [2], 1), EPSILON);
        Assert.AreEqual(0, Fitness.Zygosity(testNoEffect, [2], 0), EPSILON);

        Gene[] testMissing = [MakeGene("chr1", 0.1, 0)];
        Assert.AreEqual(0, Fitness.Zygosity(testMissing, [0], 1), EPSILON);
        Assert.AreEqual(1, Fitness.Zygosity(testMissing, [0],0), EPSILON);

        Gene[] testHaplosufficient = [MakeGene("chr1", 0.1, 0)];
        Assert.AreEqual(1, Fitness.Zygosity(testHaplosufficient, [1],1), EPSILON);
        Assert.AreEqual(0, Fitness.Zygosity(testHaplosufficient, [1],0), EPSILON);
        
        Gene[] testList = [
            MakeGene("chr1", 0.1, 0), 
            MakeGene("chr2", 0.2, 1)
        ];
        Assert.AreEqual(1, Fitness.Zygosity(testList, [1, 0], 1), EPSILON);
        Assert.AreEqual(1, Fitness.Zygosity(testList, [1, 0], 0), EPSILON);
    }

    [Test]
    public void TestEmptyTsgOgAut([Values] SexType sex, [Values(0,1)] int refId)
    {
        
        var geneList = _refs[refId].SexGeneLists[(int) sex];
        var emptyGenes = _refs[refId].GetInitialGeneCounts(sex, true);
        foreach (int geneLT in Enumerable.Range(0, geneList.Count))
        {
            Assert.AreEqual(Fitness.TsgOgTerm(geneList[geneLT], emptyGenes[geneLT]), 0, EPSILON);
        }
    }
    
    [Test]
    public void TestTsgOgAut([Values] SexType sex, [Values(0,1)] int refId)
    {
        Gene[] testNoEffect = [MakeGene("chr1", 0, 0)];
        Assert.AreEqual(Fitness.TsgOgTerm(testNoEffect, [0]), 0, EPSILON);

        Gene[] testOg = [MakeGene("chr1", 0.1, 0)];
        Assert.Greater(Fitness.TsgOgTerm(testOg, [1]), 0);

        Gene[] testTsg = [MakeGene("chr1", -0.1, 0)];
        Assert.Less(Fitness.TsgOgTerm(testTsg, [1]), 0);

        Gene[] testList =
        [
            MakeGene("chr1", 0.1, 0),
            MakeGene("chr1", -0.1, 1)
        ];
        Assert.AreEqual(Fitness.TsgOgTerm(testList, [2, 2]), 0, EPSILON);
    }

    
    [Test]
    public void TestStressTerm([Values] SexType sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var kar = new Karyotype(genRef, sex);
        long genLen = genRef.SexGenomeLen[(int)sex];
        Assert.AreEqual(0, Fitness.StressTerm(genLen, kar.GenomeLen()), EPSILON);
        kar.ApplyWGD(); // Double all
        Assert.AreEqual(-1, Fitness.StressTerm(genLen, kar.GenomeLen()), EPSILON);
        kar.ApplyWGD(); // Double all
        Assert.AreEqual(-3, Fitness.StressTerm(genLen, kar.GenomeLen()), EPSILON);
    }

    [Test]
    public void TestAutosomeStressTerm([Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var karA = new Karyotype(genRef, SexType.Any);
        var karB = new Karyotype(genRef, SexType.Any);
        Assert.AreEqual(0, Fitness.StressTerm(genRef.SexGenomeLen[(int) SexType.Any], karA.GenomeLen()), EPSILON);
        karA.ApplyWGD(); // Double all
        Assert.AreEqual(-1, Fitness.StressTerm(genRef.SexGenomeLen[(int) SexType.Any], karA.GenomeLen()), EPSILON);
        karA.ApplyWGD(); // Double all
        Assert.AreEqual(-3, Fitness.StressTerm(genRef.SexGenomeLen[(int) SexType.Any], karA.GenomeLen()), EPSILON);
        foreach (int i in Enumerable.Range(0, genRef.SexGenome[(int) SexType.Any].Count / 2)) { karB.ApplyContigDeletion(i); }
        Assert.AreEqual(0, Fitness.StressTerm(genRef.SexGenomeLen[(int) SexType.Male], karB.GenomeLen()), EPSILON);
    }

    [Test]
    public void TestCalculate([Values] SexType sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var karyotype = new Karyotype(genRef, sex);
        var fit = new FitParams(0.001, 0.001, 0.00_1);
        Assert.AreEqual(1, Fitness.Calculate(karyotype, genRef, fit), EPSILON);
    }

    [Test]
    public void TestAutosomeCalculate([Values] SexType sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var karyotype = new Karyotype(genRef, sex);
        karyotype.MergeRegions();
        var fit = new FitParams(0.001, 0.001, 0.001);
        Assert.AreEqual(1, Fitness.Calculate(karyotype, genRef, fit), EPSILON);
    }

    [Test]
    public void TestReferenceFitness([Values] SexType sex, [Values(0,1)] int refId)
    {
        var kar = new Karyotype(_refs[refId], sex);
        var genes = _refs[refId].SexGeneLists[(int) sex];
        double tsg = Fitness.TsgOgTerm(genes[(int) GeneLT.TSG], kar.GeneCounts[(int) GeneLT.TSG]);
        double og = Fitness.TsgOgTerm(genes[(int) GeneLT.TSG], kar.GeneCounts[(int) GeneLT.OG]);
        Assert.AreEqual(tsg, og, EPSILON);
    }
}
