// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SimChA.Computation;
using SimChA.Data;
using SimChA.EventData;
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
        _refs = new List<GenRef> { FileIO.ReadGenRef(TestParsing.HG_19_PATH), FileIO.ReadGenRef(TestParsing.HG_38_PATH) };
    }
    
    private static Gene MakeGene(string chrNo, double deltaFitness)
        => new($"G{chrNo}", new Region(0, 50, chrNo, false), deltaFitness);

    [Test]
    public void TestEssTerm([Values] SexType sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        bool normGenes = false;

        Assert.AreEqual(0, Fitness.EssTerm(genRef, new List<(Gene, int)>(), sex, normGenes), EPSILON);
        
        var testNoEffect = new List<(Gene, int)> { (MakeGene("chr1", 0), 0) };
        Assert.AreEqual(0, Fitness.EssTerm(genRef, testNoEffect, sex, normGenes), EPSILON);
        
        var testMissing = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 0) };
        Assert.AreEqual(-0.1, Fitness.EssTerm(genRef, testMissing, sex, normGenes), EPSILON);

        var testHaplosufficient = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 1) };
        Assert.AreEqual(0, Fitness.EssTerm(genRef, testHaplosufficient, sex, normGenes), EPSILON);
        
        var testList = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 0), (MakeGene("chr2", 0.2), 0) };
        Assert.AreEqual(-0.1 + -0.2, Fitness.EssTerm(genRef, testList, sex, normGenes), EPSILON);

        var testSexChromosome = new List<(Gene, int)> { (MakeGene("chrX", 0.8), 0), (MakeGene("chrY", 0.3), 0)};
        var expectedVal = sex switch {
            SexType.Female => -0.8,
            SexType.Male => -1.1,
            _ => 0
        };
        Assert.AreEqual(expectedVal, Fitness.EssTerm(genRef, testSexChromosome, sex, normGenes), EPSILON);
    }
    

    [Test]
    public void TestZygosity()
    {
        var genRef = _refs[0];
        // Hemizygous
        Assert.AreEqual(0, Fitness.Zygosity(genRef, new List<(Gene, int)>(), 1), EPSILON);
        // Nullizygous
        Assert.AreEqual(0, Fitness.Zygosity(genRef, new List<(Gene, int)>(), 0), EPSILON);
        
        var testNoEffect = new List<(Gene, int)> { (MakeGene("chr1", 0), 2) };
        Assert.AreEqual(0, Fitness.Zygosity(genRef, testNoEffect, 1), EPSILON);
        Assert.AreEqual(0, Fitness.Zygosity(genRef, testNoEffect, 0, true), EPSILON);

        var testMissing = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 0) };
        Assert.AreEqual(0, Fitness.Zygosity(genRef, testMissing, 1), EPSILON);
        Assert.AreEqual(1, Fitness.Zygosity(genRef, testMissing, 0, true), EPSILON);

        var testHaplosufficient = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 1) };
        Assert.AreEqual(1, Fitness.Zygosity(genRef, testHaplosufficient, 1), EPSILON);
        Assert.AreEqual(0, Fitness.Zygosity(genRef, testHaplosufficient, 0, true), EPSILON);
        
        var testList = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 1), (MakeGene("chr2", 0.2), 0) };
        Assert.AreEqual(0.5, Fitness.Zygosity(genRef, testList, 1, true), EPSILON);
        Assert.AreEqual(1, Fitness.Zygosity(genRef, testList, 0), EPSILON);
    }


    [Test]
    public void TestTsgOgAut([Values] SexType sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        Assert.AreEqual(Fitness.TsgOgTerm(genRef, new List<(Gene, int)>(), sex), 0, EPSILON);
        
        var testNoEffect = new List<(Gene, int)> {(MakeGene("chr1", 0), 0)};
        Assert.AreEqual(Fitness.TsgOgTerm(genRef, testNoEffect, sex), 0, EPSILON);

        var testOg = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 1)};
        Assert.Greater(Fitness.TsgOgTerm(genRef, testOg, sex), 0);

        var testTsg = new List<(Gene, int)> {(MakeGene("chr1", -0.1), 1)};
        Assert.Less(Fitness.TsgOgTerm(genRef, testTsg, sex), 0);

        var testList = new List<(Gene, int)> {
            (MakeGene("chr1", 0.1), 2), 
            (MakeGene("chr1", -0.1), 2), 
        };
        Assert.AreEqual(Fitness.TsgOgTerm(genRef, testList, sex), 0, EPSILON);
    }
    
    [Test]
    public void TestTsgOgSex([Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var xxGene = new List<(Gene, int)> {
            (MakeGene("chrX", 0.1), 2), 
        };
        
        var xyGeneX = new List<(Gene, int)> {
            (MakeGene("chrX", 0.1), 1), 
        };
        var xyGeneY = new List<(Gene, int)> {
            (MakeGene("chrY", 0.1), 1), 
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
    public void TestCNCalulation([Values] SexType sex,[Values(0,1)] int refId)
    {
        // Seed 14 to get chr1 delete
        var genRef = _refs[refId];
        var rnd = new Random(14);
        var karyotype = new Karyotype(genRef, sex);
        var deletion = new CNEventPars(CNEventType.ChromDeletion, 1);
        var dict = genRef.AllChrs.ToDictionary(c => c, _ => new List<Gene>());
        dict["chr1"] = new List<Gene> { MakeGene("chr1", 0.01) };
        Assert.AreEqual(Fitness.CalcCNs(dict, karyotype).FirstOrDefault(), (dict["chr1"].FirstOrDefault(), 2));
        TestKaryotype.ApplyRandomEvent(rnd, karyotype, deletion);
        Assert.AreEqual(Fitness.CalcCNs(dict, karyotype).FirstOrDefault(), (dict["chr1"].FirstOrDefault(), 1));
        TestKaryotype.ApplyRandomEvent(rnd, karyotype, deletion);
        Assert.AreEqual(Fitness.CalcCNs(dict, karyotype).FirstOrDefault(), (dict["chr1"].FirstOrDefault(), 1));
    }

    [Test]
    public void TestCalculate([Values] SexType sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var karyotype = new Karyotype(genRef, sex);
        var fit = new FitParams(0.001, 0.01, 0.000_1, true);
        Assert.AreEqual(1, Fitness.Calculate(karyotype, genRef, fit), EPSILON);
    }

    [Test]
    public void TestAutosomeCalculate([Values] SexType sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var karyotype = new Karyotype(genRef, sex);
        karyotype.MergeRegions();
        var fit = new FitParams(0.001, 0.01, 0.000_1, true);
        Assert.AreEqual(1, Fitness.Calculate(karyotype, genRef, fit), EPSILON);
    }

    [Test]
    public void TestReferenceFitness([Values] SexType sex, [Values(0, 1)] int refId, [Values(-1, 0, 1)] int myInt)
    {
        var genRef = _refs[refId];
        var karyotype = new Karyotype(genRef, sex);
        var tsgCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.TumorSuppressor], karyotype);
        var ogsCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.Oncogene], karyotype);
        double tsg = Fitness.TsgOgTerm(genRef, tsgCNs, sex);
        double og = Fitness.TsgOgTerm(genRef, ogsCNs, sex);
        Assert.AreEqual(tsg, og, EPSILON);
    }
}
