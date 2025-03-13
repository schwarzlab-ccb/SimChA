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
    public const double EPSILON = 0.00001;

    private List<GenRef> _refs;

    [SetUp]
    public void Setup()
    {
        _refs = [FileIO.ReadGenRef(TestParsing.HG_19_PATH), FileIO.ReadGenRef(TestParsing.HG_38_PATH)];
    }
    
    private static Gene MakeGene(string chrNo, double deltaFitness)
        => new(0, 50, chrNo, $"G{chrNo}", deltaFitness, GeneLT.OG);

    [Test]
    public void TestEssTerm([Values] SexType sex, [Values(0,1)] int refId)
    {
        Assert.AreEqual(0, Fitness.EssTerm( []), EPSILON);
        
        var testNoEffect = new Dictionary<Gene, int> { [MakeGene("chr1", 0)] = 0 };
        Assert.AreEqual(0, Fitness.EssTerm(testNoEffect), EPSILON);
        
        var testMissing = new Dictionary<Gene, int> { [MakeGene("chr1", 0.1)] = 0 };
        Assert.AreEqual(-0.1, Fitness.EssTerm(testMissing), EPSILON);

        var testHaplosufficient = new Dictionary<Gene, int> { [MakeGene("chr1", 0.1)] = 1 };
        Assert.AreEqual(0, Fitness.EssTerm(testHaplosufficient), EPSILON);
        
        var testList = new Dictionary<Gene, int>
        {
            [MakeGene("chr1", 0.1)] = 0,
            [MakeGene("chr2", 0.2)] = 0
        };
        Assert.AreEqual(-0.1 + -0.2, Fitness.EssTerm(testList), EPSILON);
    }

    [Test]
    public void TestZygosity()
    {
        // Hemizygous
        Assert.AreEqual(0, Fitness.Zygosity([], 1), EPSILON);
        // Nullizygous
        Assert.AreEqual(0, Fitness.Zygosity([], 0), EPSILON);
        
        var testNoEffect = new Dictionary<Gene, int> { [MakeGene("chr1", 0)] = 2 };
        Assert.AreEqual(0, Fitness.Zygosity(testNoEffect, 1), EPSILON);
        Assert.AreEqual(0, Fitness.Zygosity(testNoEffect, 0, true), EPSILON);

        var testMissing = new Dictionary<Gene, int>{ [MakeGene("chr1", 0.1)] = 0 };
        Assert.AreEqual(0, Fitness.Zygosity(testMissing, 1), EPSILON);
        Assert.AreEqual(1, Fitness.Zygosity(testMissing, 0, true), EPSILON);

        var testHaplosufficient =new Dictionary<Gene, int>{ [MakeGene("chr1", 0.1)] =  1 };
        Assert.AreEqual(1, Fitness.Zygosity(testHaplosufficient, 1), EPSILON);
        Assert.AreEqual(0, Fitness.Zygosity(testHaplosufficient, 0, true), EPSILON);
        
        var testList = new Dictionary<Gene, int>
        { 
            [MakeGene("chr1", 0.1)] = 1, 
            [MakeGene("chr2", 0.2)] = 0 
        };
        Assert.AreEqual(0.5, Fitness.Zygosity(testList, 1, true), EPSILON);
        Assert.AreEqual(1, Fitness.Zygosity(testList, 0), EPSILON);

        Assert.AreEqual(0.5, Fitness.Zygosity(testList, 1, true), EPSILON);
        Assert.AreEqual(1, Fitness.Zygosity(testList, 0), EPSILON);
    }
    
    [Test]
    public void TestTsgOgAut([Values] SexType sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        Assert.AreEqual(Fitness.TsgOgTerm(genRef, [], sex), 0, EPSILON);
        
        var testNoEffect = new Dictionary<Gene, int>{[MakeGene("chr1", 0)] = 0};
        Assert.AreEqual(Fitness.TsgOgTerm(genRef, testNoEffect, sex), 0, EPSILON);

        var testOg = new Dictionary<Gene, int>{[MakeGene("chr1", 0.1)] = 1};
        Assert.Greater(Fitness.TsgOgTerm(genRef, testOg, sex), 0);

        var testTsg = new Dictionary<Gene, int>{[MakeGene("chr1", -0.1)] = 1};
        Assert.Less(Fitness.TsgOgTerm(genRef, testTsg, sex), 0);

        var testList = new Dictionary<Gene, int>
        {
            [MakeGene("chr1", 0.1)] = 2, 
            [MakeGene("chr1", -0.1)] = 2 
        };
        Assert.AreEqual(Fitness.TsgOgTerm(genRef, testList, sex), 0, EPSILON);
    }
    
    [Test]
    public void TestTsgOgSex([Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var xxGene = new Dictionary<Gene, int> 
        {
            [MakeGene("chrX", 0.1)] = 2, 
        };
        var xyGeneX = new Dictionary<Gene, int> {
            [MakeGene("chrX", 0.1)] = 1, 
        };
        var xyGeneY = new Dictionary<Gene, int> {
            [MakeGene("chrY", 0.1)] = 1, 
        };
        
        Assert.AreEqual(
            Fitness.TsgOgTerm(genRef, xxGene, SexType.Female),
            Fitness.TsgOgTerm(genRef, xyGeneY, SexType.Male),
            EPSILON);
        
        Assert.AreEqual(
            Fitness.TsgOgTerm(genRef, xyGeneX, SexType.Male),
            Fitness.TsgOgTerm(genRef, xyGeneY, SexType.Male),
            EPSILON);
    }
    
    [Test]
    public void TestStressTerm([Values] SexType sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var kar = new Karyotype(genRef, sex);
        Assert.AreEqual(0, Fitness.StressTerm(genRef.GetGenomeLen(sex), kar.GenomeLen()), EPSILON);
        kar.ApplyWGD(); // Double all
        Assert.AreEqual(-1, Fitness.StressTerm(genRef.GetGenomeLen(sex), kar.GenomeLen()), EPSILON);
        kar.ApplyWGD(); // Double all
        Assert.AreEqual(-3, Fitness.StressTerm(genRef.GetGenomeLen(sex), kar.GenomeLen()), EPSILON);
    }

    [Test]
    public void TestAutosomeStressTerm([Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var karA = new Karyotype(genRef, SexType.Any);
        var karB = new Karyotype(genRef, SexType.Any);
        Assert.AreEqual(0, Fitness.StressTerm(genRef.GetGenomeLen(SexType.Any), karA.GenomeLen()), EPSILON);
        karA.ApplyWGD(); // Double all
        Assert.AreEqual(-1, Fitness.StressTerm(genRef.GetGenomeLen(SexType.Any), karA.GenomeLen()), EPSILON);
        karA.ApplyWGD(); // Double all
        Assert.AreEqual(-3, Fitness.StressTerm(genRef.GetGenomeLen(SexType.Any), karA.GenomeLen()), EPSILON);
        foreach (int i in Enumerable.Range(0, genRef.ChrCount(SexType.Any, false))) { karB.ApplyContigDeletion(i); }
        Assert.AreEqual(0, Fitness.StressTerm(genRef.GetGenomeLen(SexType.Male), karB.GenomeLen()), EPSILON);
    }
    

    [Test]
    public void TestCalculate([Values] SexType sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var karyotype = new Karyotype(genRef, sex);
        var fit = new FitParams(0.001, 0.01, 0.000_1, true, true);
        Assert.AreEqual(1, Fitness.Calculate(karyotype, genRef, fit), EPSILON);
    }

    [Test]
    public void TestAutosomeCalculate([Values] SexType sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var karyotype = new Karyotype(genRef, sex);
        karyotype.MergeRegions();
        var fit = new FitParams(0.001, 0.01, 0.000_1, true, true);
        Assert.AreEqual(1, Fitness.Calculate(karyotype, genRef, fit), EPSILON);
    }

    [Test]
    public void TestReferenceFitness([Values] SexType sex, [Values(0, 1)] int refId, [Values(-1, 0, 1)] int myInt)
    {
        var genRef = _refs[refId];
        var kar = new Karyotype(genRef, sex);
        double tsg = Fitness.TsgOgTerm(genRef, kar.GeneCounts[GeneLT.TSG], sex);
        double og = Fitness.TsgOgTerm(genRef, kar.GeneCounts[GeneLT.OG], sex);
        Assert.AreEqual(tsg, og, EPSILON);
    }
}
