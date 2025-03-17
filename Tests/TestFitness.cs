// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

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

    private List<GenRef> _refs;

    [SetUp]
    public void Setup()
    {
        _refs = [FileIO.ReadGenRef(TestParsing.HG_19_PATH), FileIO.ReadGenRef(TestParsing.HG_38_PATH)];
    }
    
    private static Gene MakeGene(string chrNo, double deltaFitness, int index)
        => new(0, 50, chrNo, GeneLT.OG, index, deltaFitness);

    [Test]
    public void TestEssTerm([Values] SexType sex, [Values(0,1)] int refId)
    {
        Assert.AreEqual(0, Fitness.EssTerm( [], []), EPSILON);
        
        List<Gene> testNoEffect = [MakeGene("chr1", 0, 0)];
        Assert.AreEqual(0, Fitness.EssTerm(testNoEffect, [0]), EPSILON);
        
        List<Gene> testMissing = [MakeGene("chr1", 0.1, 0)];
        Assert.AreEqual(-0.1, Fitness.EssTerm(testMissing, [0]), EPSILON);

        List<Gene> testHaplosufficient = [MakeGene("chr1", 0.1, 0)];
        Assert.AreEqual(0, Fitness.EssTerm(testHaplosufficient, [1]), EPSILON);

        List<Gene> testList = [
            MakeGene("chr1", 0.1, 0), 
            MakeGene("chr2", 0.2, 1)
        ];
        Assert.AreEqual(-0.1 + -0.2, Fitness.EssTerm(testList, [0, 0]), EPSILON);
    }

    [Test]
    public void TestZygosity()
    {
        // Hemizygous
        Assert.AreEqual(0, Fitness.Zygosity([], [], 1), EPSILON);
        // Nullizygous
        Assert.AreEqual(0, Fitness.Zygosity([], [], 0), EPSILON);
        
        List<Gene> testNoEffect = [MakeGene("chr1", 0, 0)];
        Assert.AreEqual(0, Fitness.Zygosity(testNoEffect, [2], 1), EPSILON);
        Assert.AreEqual(0, Fitness.Zygosity(testNoEffect, [2], 0, true), EPSILON);

        List<Gene> testMissing = [MakeGene("chr1", 0.1, 0)];
        Assert.AreEqual(0, Fitness.Zygosity(testMissing, [0], 1), EPSILON);
        Assert.AreEqual(1, Fitness.Zygosity(testMissing, [0],0, true), EPSILON);

        List<Gene> testHaplosufficient = [MakeGene("chr1", 0.1, 0)];
        Assert.AreEqual(1, Fitness.Zygosity(testHaplosufficient, [1],1), EPSILON);
        Assert.AreEqual(0, Fitness.Zygosity(testHaplosufficient, [1],0, true), EPSILON);
        
        List<Gene> testList = [
            MakeGene("chr1", 0.1, 0), 
            MakeGene("chr2", 0.2, 1)
        ];
        Assert.AreEqual(0.5, Fitness.Zygosity(testList, [1, 0], 1, true), EPSILON);
        Assert.AreEqual(1, Fitness.Zygosity(testList, [1, 0], 0), EPSILON);
    }

    [Test]
    public void TestEmptyTsgOgAut([Values] SexType sex, [Values(0,1)] int refId)
    {
        
        var geneList = _refs[refId].GeneLists[(int) sex];
        var emptyGenes = _refs[refId].GetInitialGeneCounts(sex, true);
        for (int geneLT = 0; geneLT < 3; geneLT++)
        {
            Assert.AreEqual(Fitness.TsgOgTerm(geneList[geneLT], emptyGenes[geneLT]), 0, EPSILON);
        }
    }
    
    [Test]
    public void TestTsgOgAut([Values] SexType sex, [Values(0,1)] int refId)
    {
        List<Gene> testNoEffect = [MakeGene("chr1", 0, 0)];
        Assert.AreEqual(Fitness.TsgOgTerm(testNoEffect, [0]), 0, EPSILON);

        List<Gene> testOg = [MakeGene("chr1", 0.1, 0)];
        Assert.Greater(Fitness.TsgOgTerm(testOg, [1]), 0);

        List<Gene> testTsg = [MakeGene("chr1", -0.1, 0)];
        Assert.Less(Fitness.TsgOgTerm(testTsg, [1]), 0);

        List<Gene> testList =
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
        long genLen = genRef.GenomeLens[(int)sex];
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
        Assert.AreEqual(0, Fitness.StressTerm(genRef.GenomeLens[(int) SexType.Any], karA.GenomeLen()), EPSILON);
        karA.ApplyWGD(); // Double all
        Assert.AreEqual(-1, Fitness.StressTerm(genRef.GenomeLens[(int) SexType.Any], karA.GenomeLen()), EPSILON);
        karA.ApplyWGD(); // Double all
        Assert.AreEqual(-3, Fitness.StressTerm(genRef.GenomeLens[(int) SexType.Any], karA.GenomeLen()), EPSILON);
        foreach (int i in Enumerable.Range(0, genRef.Genomes[(int) SexType.Any].Count / 2)) { karB.ApplyContigDeletion(i); }
        Assert.AreEqual(0, Fitness.StressTerm(genRef.GenomeLens[(int) SexType.Male], karB.GenomeLen()), EPSILON);
    }

    [Test]
    public void TestCalculate([Values] SexType sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var karyotype = new Karyotype(genRef, sex);
        var fit = new FitParams(0.001, 0.001, 0.00_1, true, true);
        Assert.AreEqual(1, Fitness.Calculate(karyotype, genRef, fit), EPSILON);
    }

    [Test]
    public void TestAutosomeCalculate([Values] SexType sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var karyotype = new Karyotype(genRef, sex);
        karyotype.MergeRegions();
        var fit = new FitParams(0.001, 0.001, 0.001, true, true);
        Assert.AreEqual(1, Fitness.Calculate(karyotype, genRef, fit), EPSILON);
    }

    [Test]
    public void TestReferenceFitness([Values] SexType sex, [Values(0,1)] int refId)
    {
        var kar = new Karyotype(_refs[refId], sex);
        var genes = _refs[refId].GeneLists[(int) sex];
        double tsg = Fitness.TsgOgTerm(genes[(int) GeneLT.TSG], kar.GeneCounts[(int) GeneLT.TSG]);
        double og = Fitness.TsgOgTerm(genes[(int) GeneLT.TSG], kar.GeneCounts[(int) GeneLT.OG]);
        Assert.AreEqual(tsg, og, EPSILON);
    }
}
