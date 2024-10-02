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
    public void TestEssTerm([Values(0,1)] int refId)
    {
        var genRef = _refs[refId];

        Assert.AreEqual(0, Fitness.EssTerm(genRef, new List<(Gene, int)>(), SexEnum.Female));
        
        var testNoEffect = new List<(Gene, int)> { (MakeGene("chr1", 0), 0) };
        Assert.AreEqual(0, Fitness.EssTerm(genRef, testNoEffect, SexEnum.Female));
        
        var testMissing = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 0) };
        Assert.AreEqual(-0.1, Fitness.EssTerm(genRef, testMissing, SexEnum.Female));

        var testHaplosufficient = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 1) };
        Assert.AreEqual(0, Fitness.EssTerm(genRef, testHaplosufficient, SexEnum.Female));
        
        var testList = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 0), (MakeGene("chr2", 0.2), 0) };
        Assert.AreEqual(-0.1 + -0.2, Fitness.EssTerm(genRef, testList, SexEnum.Female));

        var testSexChromosome = new List<(Gene, int)> { (MakeGene("chrX", 0.8), 0), (MakeGene("chrY", 0.3), 0)};
        Assert.AreEqual(0.0, Fitness.EssTerm(genRef, testSexChromosome, SexEnum.None));

        Assert.AreEqual(-0.8, Fitness.EssTerm(genRef, testSexChromosome, SexEnum.Female));
        Assert.AreEqual(-1.1, Fitness.EssTerm(genRef, testSexChromosome, SexEnum.Male));
    }

    [Test]
    public void TestEssTermHaploinsufficiency([Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        genRef.IncludeSexChromosomes = false;
        Assert.AreEqual(0, Fitness.EssTerm(genRef, new List<(Gene, int)>(), SexEnum.Female, false, true));
        
        var testNoEffect = new List<(Gene, int)> { (MakeGene("chr1", 0), 0) };
        Assert.AreEqual(0, Fitness.EssTerm(genRef, testNoEffect, SexEnum.Female, false, true));
        
        var testMissing = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 0) };
        Assert.AreEqual(-0.2, Fitness.EssTerm(genRef, testMissing, SexEnum.Female, false, true));

        var testHaploinsufficient = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 1) };
        Assert.AreEqual(-0.1, Fitness.EssTerm(genRef, testHaploinsufficient, SexEnum.Female, false, true));
        
        var testList = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 0), (MakeGene("chr2", 0.2), 0) };
        Assert.AreEqual(-0.2 + -0.4, Fitness.EssTerm(genRef, testList, SexEnum.Female, false, true));

        var testSexChromosome = new List<(Gene, int)> { (MakeGene("chrX", 0.5), 0), (MakeGene("chrY", 0.5), 0)};
        Assert.AreEqual(0.0 + 0.0, Fitness.EssTerm(genRef, testSexChromosome, SexEnum.Female, false, true));

        genRef.IncludeSexChromosomes = true;
        Assert.AreEqual(-1.0 +  0.0, Fitness.EssTerm(genRef, testSexChromosome, SexEnum.Female, false, true));
        Assert.AreEqual(-0.5 + -0.5, Fitness.EssTerm(genRef, testSexChromosome, SexEnum.Male, false, true));
    }

    [Test]
    public void TestTsgOgTerm([Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        Assert.AreEqual(0, Fitness.TsgOgTerm(genRef, new List<(Gene, int)>(), SexEnum.Female));
        
        var testNoEffect = new List<(Gene, int)> {(MakeGene("chr1", 0), 0)};
        Assert.AreEqual(0, Fitness.TsgOgTerm(genRef, testNoEffect, SexEnum.Female));

        var testMissing = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 1)};
        Assert.AreEqual(-0.1, Fitness.TsgOgTerm(genRef, testMissing, SexEnum.Female));

        var testMissingTwice = new List<(Gene, int)> {(MakeGene("chr1", 0.1), 0)};
        Assert.AreEqual(-0.2, Fitness.TsgOgTerm(genRef, testMissingTwice, SexEnum.Female));

        var testList = new List<(Gene, int)> {
            (MakeGene("chr1", 0.1), 1), 
            (MakeGene("chr1", 0.2), 0), 
            (MakeGene("chr1", 0.3), 2)
        };
        Assert.AreEqual(-0.2 - 0.2 - 0.1, Fitness.TsgOgTerm(genRef, testList, SexEnum.Female));

        testList = new List<(Gene, int)>{
            (MakeGene("chrX", 0.1), 2)
        };
        Assert.AreEqual(0, Fitness.TsgOgTerm(genRef, testList, SexEnum.Female));
        Assert.AreEqual(0.1, Fitness.TsgOgTerm(genRef, testList, SexEnum.Male));
        testList = new List<(Gene, int)>{
            (MakeGene("chrX", 0.2), 2),
            (MakeGene("chrY", 0.1), 2)
        };
        Assert.AreEqual(0.2, Fitness.TsgOgTerm(genRef, testList, SexEnum.Female));
        Assert.AreEqual(0.2 + 0.1, Fitness.TsgOgTerm(genRef, testList, SexEnum.Male));
        // Autosomes only
        genRef.IncludeSexChromosomes = false;
                testList = new List<(Gene, int)>{
            (MakeGene("chrX", 0.1), 2)
        };
        Assert.AreEqual(0, Fitness.TsgOgTerm(genRef, testList, SexEnum.Female));
        Assert.AreEqual(0, Fitness.TsgOgTerm(genRef, testList, SexEnum.Male));
        testList = new List<(Gene, int)>{
            (MakeGene("chrX", 0.2), 2),
            (MakeGene("chrY", 0.1), 2)
        };
        Assert.AreEqual(0, Fitness.TsgOgTerm(genRef, testList, SexEnum.Female));
        Assert.AreEqual(0, Fitness.TsgOgTerm(genRef, testList, SexEnum.Male));

    }
    
    [Test]
    public void TestStressTerm([Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var xxKaryotype = new Karyotype(genRef, SexEnum.Female);
        var xyKaryotype = new Karyotype(genRef, SexEnum.Male);
        Assert.AreEqual(0, Fitness.StressTerm(genRef.GetGenomeLen(SexEnum.Female), xxKaryotype.GenomeLen()), EPSILON);
        Assert.AreEqual(0, Fitness.StressTerm(genRef.GetGenomeLen(SexEnum.Male), xyKaryotype.GenomeLen()), EPSILON);
        Assert.AreNotEqual(0, Fitness.StressTerm(genRef.GetGenomeLen(SexEnum.Male), xxKaryotype.GenomeLen()));
        xxKaryotype.ApplyWGD(); // Double all
        Assert.AreEqual(-1, Fitness.StressTerm(genRef.GetGenomeLen(SexEnum.Female), xxKaryotype.GenomeLen()), EPSILON);
        xxKaryotype.ApplyWGD(); // Double all
        Assert.AreEqual(-3, Fitness.StressTerm(genRef.GetGenomeLen(SexEnum.Female), xxKaryotype.GenomeLen()), EPSILON);
        foreach (int i in Enumerable.Range(0, genRef.ChrCount(SexEnum.None, false))) { xyKaryotype.ApplyContigDeletion(i); }
        Assert.AreEqual(0, Fitness.StressTerm(genRef.GetGenomeLen(SexEnum.Male), xyKaryotype.GenomeLen()), EPSILON);
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
    public void TestCNCalulation([Values(0,1)] int refId)
    {
        // Seed 14 to get chr1 delete
        var genRef = _refs[refId];
        var rnd = new Random(14);
        var karyotype = new Karyotype(genRef, SexEnum.Female);
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
    public void TestCalculate([Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var karyotype = new Karyotype(genRef, SexEnum.Female);
        var fit = new FitnessParams(0.001f, 0.01f, 0.000_1f, 1f);
        Assert.AreEqual(1, Fitness.Calculate(karyotype, genRef, fit), EPSILON);
        // TODO: Test the linear combination
    }

    [Test]
    public void TestCalculateFromComponents([Values] SexEnum sex, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var karyotype = new Karyotype(genRef, sex);
        var fit = new FitnessParams(0.001f, 0.01f, 0.000_1f, 1f);
        double stress = Fitness.StressTerm(genRef.GetGenomeLen(sex), karyotype.GenomeLen());
        double tsg = -Fitness.TsgOgTerm(genRef, Fitness.CalcCNs(genRef.GeneLists[GeneListType.TumorSuppressor], karyotype), sex);
        double og = Fitness.TsgOgTerm(genRef, Fitness.CalcCNs(genRef.GeneLists[GeneListType.Oncogene], karyotype), sex);
        double ess = Fitness.EssTerm(genRef, Fitness.CalcCNs(genRef.GeneLists[GeneListType.Essentiality], karyotype), sex);
        double total = 1 + (stress*fit.Stress + (tsg + og)*fit.TsgOg + ess*fit.Essentiality) * fit.TotalStrength;
        Assert.AreEqual(total, Fitness.CalculateFromComponents(stress, tsg+og, ess, fit), EPSILON);
    }

    [Test]
    public void TestAutosomeCalculate([Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        genRef.IncludeSexChromosomes = false;
        var karyotype = new Karyotype(genRef, SexEnum.Female);
        karyotype.MergeRegions();
        var fit = new FitnessParams(0.001f, 0.01f, 0.000_1f, 1f);
        Assert.AreEqual(1, Fitness.Calculate(karyotype, genRef, fit), EPSILON);
        // TODO: Test the linear combination
    }
    
    [Test]
    public void TestReferenceFitness([Values] SexEnum sexXX, [Values(0,1)] int refId, [Values(-1, 0, 1)] int myInt)
    {
        var genRef = _refs[refId];
        var karyotype = new Karyotype(genRef, sexXX);
        var tsgCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.TumorSuppressor], karyotype);
        double tsg = Fitness.TsgOgTerm(genRef, tsgCNs, sexXX);
        Assert.AreEqual(0, tsg, EPSILON);
        var ogsCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.Oncogene], karyotype);
        double og = Fitness.TsgOgTerm(genRef, ogsCNs, sexXX);;
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
