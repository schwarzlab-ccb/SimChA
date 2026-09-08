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

    // The exponent the default FitParams carries; the graded assertions below are written against it.
    private const double K = 2.0;
    // Large enough that 0.5^k underflows the tolerance, i.e. the CN==0-only step EssTerm replaced.
    private const double K_STEP = 100.0;

    [Test]
    public void TestEssTerm([Values] SexType sex, [Values(0,1)] int refId)
    {
        Assert.AreEqual(0, Fitness.EssTerm([], [], K), EPSILON);

        // A zero-scored gene contributes nothing whatever its copy number.
        Gene[] testNoEffect = [MakeGene("chr1", 0, 0)];
        Assert.AreEqual(0, Fitness.EssTerm(testNoEffect, [0], K), EPSILON);

        // Homozygous loss costs the gene's whole score, independently of the exponent.
        Gene[] testMissing = [MakeGene("chr1", 0.1, 0)];
        Assert.AreEqual(-0.1, Fitness.EssTerm(testMissing, [0], K), EPSILON);
        Assert.AreEqual(-0.1, Fitness.EssTerm(testMissing, [0], K_STEP), EPSILON);

        // One copy left: charged 0.5^k, where the old model charged nothing at all.
        Gene[] testHemizygous = [MakeGene("chr1", 0.1, 0)];
        Assert.AreEqual(-0.1 * 0.25, Fitness.EssTerm(testHemizygous, [1], K), EPSILON);
        Assert.AreEqual(-0.1 * 0.5, Fitness.EssTerm(testHemizygous, [1], 1.0), EPSILON);
        Assert.AreEqual(0, Fitness.EssTerm(testHemizygous, [1], K_STEP), EPSILON);

        // At or above the diploid reference there is no penalty, and gains are not rewarded.
        Assert.AreEqual(0, Fitness.EssTerm(testHemizygous, [2], K), EPSILON);
        Assert.AreEqual(0, Fitness.EssTerm(testHemizygous, [5], K), EPSILON);

        Gene[] testList = [
            MakeGene("chr1", 0.1, 0), 
            MakeGene("chr2", 0.2, 1)
        ];
        Assert.AreEqual(-0.1 + -0.2, Fitness.EssTerm(testList, [0, 0], K), EPSILON);
        Assert.AreEqual(-0.1 + -0.2 * 0.25, Fitness.EssTerm(testList, [0, 1], K), EPSILON);
    }

    [Test]
    public void TestDosageLoss()
    {
        // The exponent only ever reshapes the interior; the endpoints are pinned.
        foreach (double k in new[] { 0.5, 1.0, 2.0, K_STEP })
        {
            Assert.AreEqual(1.0, Fitness.DosageLoss(0, k), EPSILON);
            Assert.AreEqual(0.0, Fitness.DosageLoss(2, k), EPSILON);
            Assert.AreEqual(0.0, Fitness.DosageLoss(9, k), EPSILON);
        }
        Assert.AreEqual(Math.Pow(0.5, 0.5), Fitness.DosageLoss(1, 0.5), EPSILON);
        Assert.AreEqual(0.5,  Fitness.DosageLoss(1, 1.0), EPSILON);
        Assert.AreEqual(0.25, Fitness.DosageLoss(1, 2.0), EPSILON);
        Assert.AreEqual(0.0,  Fitness.DosageLoss(1, K_STEP), EPSILON);

        // Monotone in the exponent at the one point that is free to move.
        Assert.Less(Fitness.DosageLoss(1, 3.0), Fitness.DosageLoss(1, 2.0));
    }

    [Test]
    public void TestAcceptProb()
    {
        // delta is the fitness gain at which a proposal is accepted half the time.
        Assert.AreEqual(0.5, Fitness.AcceptProb(0.4, 0.4), EPSILON);
        Assert.AreEqual(0.5, Fitness.AcceptProb(0.0, 0.0), EPSILON);

        // Strictly monotone in the fitness change, and strictly inside (0, 1) -- never clipped, so
        // the derivative with respect to delta is nonzero everywhere. This is the whole point of
        // the rule: Metropolis returned exactly 1 for every one of these.
        // The range spans the fitness changes actually seen in a fit: the September 2026 spice run's
        // accepted events ran from -8.8 to +20.8, and 34% of them sat above delta, where Metropolis
        // returned exactly 1 and every gradient vanished.
        double prev = 0.0;
        foreach (double dF in new[] { -8.8, -5.0, -1.0, -0.1, 0.0, 0.1, 1.0, 5.0, 20.8 })
        {
            double p = Fitness.AcceptProb(dF, 0.1);
            Assert.Greater(p, prev);
            Assert.Greater(p, 0.0);
            Assert.Less(p, 1.0);
            prev = p;
        }

        // Strictly below 1 only until the exponential underflows: 1/(1+exp(-x)) rounds to 1.0 in
        // float64 once x is past ~37. No fitted run comes near that -- the observed maximum is 20.8
        // -- but the bound is the honest limit of "never saturates", so it is pinned here.
        Assert.Less(Fitness.AcceptProb(36.0, 0.1), 1.0);
        Assert.AreEqual(1.0, Fitness.AcceptProb(40.0, 0.1), EPSILON);

        // Deep in the deleterious tail it agrees with the Metropolis rule it replaces.
        foreach (double dF in new[] { -10.0, -20.0 })
        {
            Assert.AreEqual(Math.Exp(dF - 0.1), Fitness.AcceptProb(dF, 0.1), EPSILON);
        }

        // No overflow or NaN at the extremes.
        Assert.AreEqual(0.0, Fitness.AcceptProb(-1e6, 0.1), EPSILON);
        Assert.AreEqual(1.0, Fitness.AcceptProb(1e6, 0.1), EPSILON);
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
