// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;
using SimChA.IO;
using SimChA.Simulation;

namespace Tests;

[TestFixture]
public class TestFitness
{
    public const double EPSILON = 0.0000000001;

    private List<GenRef> _refs;

    [SetUp]
    public void Setup()
    {
        _refs = new List<GenRef> { FileIO.GetGenRef(TestIO.HG_19_PATH), FileIO.GetGenRef(TestIO.HG_38_PATH) };
    }
    
    private static Gene MakeGene(string chrNo, double deltaFitness)
        => new($"G{chrNo}", new Region(0, 50, chrNo, false), deltaFitness);

    [Test]
    public void TestEssTerm([Values] SexEnum sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];

        Assert.AreEqual(0, Fitness.EssTerm(genRef, new List<(Gene, int)>(), sex), EPSILON);
        
        var testNoEffect = new List<(Gene, int)> { (MakeGene("chr1", 0), 0) };
        Assert.AreEqual(0, Fitness.EssTerm(genRef, testNoEffect, sex), EPSILON);
        
        var testMissing = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 0) };
        Assert.AreEqual(-0.1, Fitness.EssTerm(genRef, testMissing, sex), EPSILON);

        var testHaplosufficient = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 1) };
        Assert.AreEqual(0, Fitness.EssTerm(genRef, testHaplosufficient, sex), EPSILON);
        
        var testList = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 0), (MakeGene("chr2", 0.2), 0) };
        Assert.AreEqual(-0.1 + -0.2, Fitness.EssTerm(genRef, testList, sex), EPSILON);

        var testSexChromosome = new List<(Gene, int)> { (MakeGene("chrX", 0.8), 0), (MakeGene("chrY", 0.3), 0)};
        var expectedVal = sex switch {
            SexEnum.Female => -0.8,
            SexEnum.Male => -1.1,
            _ => 0
        };
        Assert.AreEqual(expectedVal, Fitness.EssTerm(genRef, testSexChromosome, sex), EPSILON);
    }

    [Test]
    public void TestEssTermHaploinsufficiency([Values] SexEnum sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        Assert.AreEqual(0, Fitness.EssTerm(genRef, new List<(Gene, int)>(), sex, false, true));
        
        var testNoEffect = new List<(Gene, int)> { (MakeGene("chr1", 0), 0) };
        Assert.AreEqual(0, Fitness.EssTerm(genRef, testNoEffect, sex, false, true), EPSILON);
        
        var testMissing = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 0) };
        Assert.AreEqual(-0.2, Fitness.EssTerm(genRef, testMissing, sex, false, true), EPSILON);

        var testHaploinsufficient = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 1) };
        Assert.AreEqual(-0.1, Fitness.EssTerm(genRef, testHaploinsufficient, sex, false, true), EPSILON);
        
        var testList = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 0), (MakeGene("chr2", 0.2), 0) };
        Assert.AreEqual(-0.2 + -0.4, Fitness.EssTerm(genRef, testList, sex, false, true), EPSILON);

        var testSexChromosome = new List<(Gene, int)> { (MakeGene("chrX", 0.7), 0), (MakeGene("chrY", 0.5), 0)};
        var expectedVal = sex switch {
            SexEnum.Female => -1.4,
            SexEnum.Male => -1.2,
            _ => 0
        };
        Assert.AreEqual(expectedVal, Fitness.EssTerm(genRef, testSexChromosome, sex, false, true), EPSILON);
    }

    [Test]
    public void TestTsgOgTerm([Values] SexEnum sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        Assert.AreEqual(0, Fitness.TsgOgTerm(genRef, new List<(Gene, int)>(), sex), EPSILON);
        
        var testNoEffect = new List<(Gene, int)> {(MakeGene("chr1", 0), 0)};
        Assert.AreEqual(0, Fitness.TsgOgTerm(genRef, testNoEffect, sex), EPSILON);

        var testMissing = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 1)};
        Assert.AreEqual(-0.1, Fitness.TsgOgTerm(genRef, testMissing, sex), EPSILON);

        var testMissingTwice = new List<(Gene, int)> {(MakeGene("chr1", 0.1), 0)};
        Assert.AreEqual(-0.2, Fitness.TsgOgTerm(genRef, testMissingTwice, sex), EPSILON);

        var testList = new List<(Gene, int)> {
            (MakeGene("chr1", 0.1), 1), 
            (MakeGene("chr1", 0.2), 0), 
            (MakeGene("chr1", 0.3), 2)
        };
        Assert.AreEqual(-0.2 - 0.2 - 0.1, Fitness.TsgOgTerm(genRef, testList, sex), EPSILON);

        testList = new List<(Gene, int)>{
            (MakeGene("chrX", 0.1), 2)
        };
        var expectedVal = sex switch {
            SexEnum.Female => 0,
            SexEnum.Male => 0.1,
            _ => 0
        };
        Assert.AreEqual(expectedVal, Fitness.TsgOgTerm(genRef, testList, sex), EPSILON);
        testList = new List<(Gene, int)>{
            (MakeGene("chrX", 0.2), 2),
            (MakeGene("chrY", 0.1), 2)
        };
        expectedVal = sex switch {
            SexEnum.Female => 0,
            SexEnum.Male => 0.3,
            _ => 0
        };
        Assert.AreEqual(expectedVal, Fitness.TsgOgTerm(genRef, testList, sex), EPSILON);
    }
    
    [Test]
    public void TestStressTerm([Values] SexEnum sex, [Values(0,1)] int refId)
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
        var karA = new Karyotype(genRef, SexEnum.None);
        var karB = new Karyotype(genRef, SexEnum.None);
        Assert.AreEqual(0, Fitness.StressTerm(genRef.GetGenomeLen(SexEnum.None), karA.GenomeLen()), EPSILON);
        karA.ApplyWGD(); // Double all
        Assert.AreEqual(-1, Fitness.StressTerm(genRef.GetGenomeLen(SexEnum.None), karA.GenomeLen()), EPSILON);
        karA.ApplyWGD(); // Double all
        Assert.AreEqual(-3, Fitness.StressTerm(genRef.GetGenomeLen(SexEnum.None), karA.GenomeLen()), EPSILON);
        foreach (int i in Enumerable.Range(0, genRef.ChrCount(SexEnum.None, false))) { karB.ApplyContigDeletion(i); }
        Assert.AreEqual(0, Fitness.StressTerm(genRef.GetGenomeLen(SexEnum.Male), karB.GenomeLen()), EPSILON);
    }
    
    [Test]
    public void TestCNCalulation([Values] SexEnum sex,[Values(0,1)] int refId)
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
    public void TestCalculate([Values] SexEnum sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var karyotype = new Karyotype(genRef, sex);
        var fit = new FitnessParams(0.001f, 0.01f, 0.000_1f);
        Assert.AreEqual(1, Fitness.Calculate(karyotype, genRef, fit), EPSILON);
        // TODO: Test the linear combination
    }

    [Test]
    public void TestCalculateFromComponents([Values] SexEnum sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var karyotype = new Karyotype(genRef, sex);
        var fit = new FitnessParams(0.001f, 0.01f, 0.000_1f);
        double stress = Fitness.StressTerm(genRef.GetGenomeLen(sex), karyotype.GenomeLen());
        double tsg = -Fitness.TsgOgTerm(genRef, Fitness.CalcCNs(genRef.GeneLists[GeneListType.TumorSuppressor], karyotype), sex);
        double og = Fitness.TsgOgTerm(genRef, Fitness.CalcCNs(genRef.GeneLists[GeneListType.Oncogene], karyotype), sex);
        double ess = Fitness.EssTerm(genRef, Fitness.CalcCNs(genRef.GeneLists[GeneListType.Essentiality], karyotype), sex);
        double total = 1 + (stress*fit.Stress + (tsg + og)*fit.TsgOg + ess*fit.Essentiality);
        Assert.AreEqual(total, Fitness.CalculateFromComponents(stress, tsg+og, ess, fit), EPSILON);
    }

    [Test]
    public void TestAutosomeCalculate([Values] SexEnum sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var karyotype = new Karyotype(genRef, sex);
        karyotype.MergeRegions();
        var fit = new FitnessParams(0.001f, 0.01f, 0.000_1f);
        Assert.AreEqual(1, Fitness.Calculate(karyotype, genRef, fit), EPSILON);
        // TODO: Test the linear combination
    }
    
    [Test]
    public void TestReferenceFitness([Values] SexEnum sex, [Values(0,1)] int refId, [Values(-1, 0, 1)] int myInt)
    {
        var genRef = _refs[refId];
        var karyotype = new Karyotype(genRef, sex);
        var tsgCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.TumorSuppressor], karyotype);
        double tsg = Fitness.TsgOgTerm(genRef, tsgCNs, sex);
        Assert.AreEqual(0, tsg, EPSILON);
        var ogsCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.Oncogene], karyotype);
        double og = Fitness.TsgOgTerm(genRef, ogsCNs, sex);;
        Assert.AreEqual(0, og, EPSILON);
    }

    [Test]
    public void TestGetPresentGenes([Values] bool useTSG, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var selectList = genRef.GeneLists[useTSG ? GeneListType.TumorSuppressor : GeneListType.Oncogene];
        var contigs = genRef.GetGenotype(SexEnum.Male).Select(region => new Contig(region)).ToList();
        foreach (string chrNo in genRef.AllChrs)
        {
            int chrToCont = chrNo != "chrY" ? genRef.AllChrs.FindIndex(c => c == chrNo) : 45;
            int contigCount = contigs[chrToCont].GetPresentGenes(selectList).Count;
            int chrCount = selectList[chrNo].Count;
            Assert.AreEqual(chrCount, contigCount);
        }
    }
    
    [Test]
    public void TestTsgOgSum([Values] SexEnum sex, [Values] bool useTSG, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var selectList = genRef.GeneLists[useTSG ? GeneListType.TumorSuppressor : GeneListType.Oncogene];
        selectList = sex switch {
            SexEnum.Female => selectList.Where(pair => pair.Key != "chrY").ToDictionary(pair => pair.Key, pair => pair.Value),
            SexEnum.Male => selectList,
            _ => selectList.Where(pair => pair.Key != "chrX" && pair.Key != "chrY").ToDictionary(pair => pair.Key, pair => pair.Value)
        };
        double sumHap1 = selectList.Where(pair => pair.Key != "chrY").Sum(pair => pair.Value.Sum(g => g.DeltaFitness));
        string missingChr = sex == SexEnum.Female ? "chrY" : "chrX";
        double sumHap2 = selectList.Where(pair => pair.Key != missingChr).Sum(pair => pair.Value.Sum(g => g.DeltaFitness));
        double total = sumHap1 + sumHap2;
        var karyotype = new Karyotype(genRef, sex);
        karyotype.ApplyWGD();
        var cnList = Fitness.CalcCNs(selectList, karyotype);
        double sum = Fitness.TsgOgTerm(genRef, cnList, sex);
        Assert.AreEqual(total, sum, EPSILON);
    }
}
