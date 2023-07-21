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
    public void TestEssTerm()
    {
        Assert.AreEqual(0, Fitness.EssTerm(new List<(Gene, int)>()));
        
        var testNoEffect = new List<(Gene, int)> { (MakeGene("chr1", 0), 0) };
        Assert.AreEqual(0, Fitness.EssTerm(testNoEffect));
        
        var testMissing = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 0) };
        Assert.AreEqual(-0.1, Fitness.EssTerm(testMissing));

        var testHaplosufficient = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 1) };
        Assert.AreEqual(0, Fitness.EssTerm(testHaplosufficient));
        
        var testList = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 0), (MakeGene("chr2", 0.2), 0) };
        Assert.AreEqual(-0.1 + -0.2, Fitness.EssTerm(testList));
    }

    [Test]
    public void TestTsgOgTerm([Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        Assert.AreEqual(0, Fitness.TsgOgTerm(genRef, new List<(Gene, int)>(), true));
        
        var testNoEffect = new List<(Gene, int)> {(MakeGene("chr1", 0), 0)};
        Assert.AreEqual(0, Fitness.TsgOgTerm(genRef, testNoEffect, true));

        var testMissing = new List<(Gene, int)> { (MakeGene("chr1", 0.1), 1)};
        Assert.AreEqual(-0.1, Fitness.TsgOgTerm(genRef, testMissing, true));

        var testMissingTwice = new List<(Gene, int)> {(MakeGene("chr1", 0.1), 0)};
        Assert.AreEqual(-0.2, Fitness.TsgOgTerm(genRef, testMissingTwice, true));

        var testList = new List<(Gene, int)> {
            (MakeGene("chr1", 0.1), 1), 
            (MakeGene("chr1", 0.2), 0), 
            (MakeGene("chr1", 0.3), 2)
        };
        Assert.AreEqual(-0.2 - 0.2 - 0.1, Fitness.TsgOgTerm(genRef, testList, true));

        testList = new List<(Gene, int)>{
            (MakeGene("chrX", 0.1), 2)
        };
        Assert.AreEqual(0, Fitness.TsgOgTerm(genRef, testList, true));
        Assert.AreEqual(0.1, Fitness.TsgOgTerm(genRef, testList, false));
        testList = new List<(Gene, int)>{
            (MakeGene("chrX", 0.2), 2),
            (MakeGene("chrY", 0.1), 2)
        };
        Assert.AreEqual(0.2, Fitness.TsgOgTerm(genRef, testList, true));
        Assert.AreEqual(0.2 + 0.1, Fitness.TsgOgTerm(genRef, testList, false));
    }
    
    [Test]
    public void TestStressTerm([Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var xxKaryotype = new Karyotype(genRef, true);
        var xyKaryotype = new Karyotype(genRef, false);
        Assert.AreEqual(0, Fitness.StressTerm(genRef.GetGenomeLen(true), xxKaryotype.GenomeLen()), EPSILON);
        Assert.AreEqual(0, Fitness.StressTerm(genRef.GetGenomeLen(false), xyKaryotype.GenomeLen()), EPSILON);
        Assert.AreNotEqual(0, Fitness.StressTerm(genRef.GetGenomeLen(false), xxKaryotype.GenomeLen()));
        xxKaryotype.ApplyWGD(); // Double all
        Assert.AreEqual(-1, Fitness.StressTerm(genRef.GetGenomeLen(true), xxKaryotype.GenomeLen()), EPSILON);
        foreach (int i in Enumerable.Range(0, genRef.ChrCount)) { xyKaryotype.ApplyContigDeletion(i); }
        Assert.AreEqual(1, Fitness.StressTerm(genRef.GetGenomeLen(true), xyKaryotype.GenomeLen()), EPSILON);
    }
    
    [Test]
    public void TestCNCalulation([Values(0,1)] int refId)
    {
        // Seed 14 to get chr1 delete
        var genRef = _refs[refId];
        var rnd = new Random(14);
        var karyotype = new Karyotype(genRef, true);
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
        var karyotype = new Karyotype(genRef, true);
        var fit = new FitnessParams(0.001f, 0.01f, 0.000_1f);
        Assert.AreEqual(1, Fitness.Calculate(karyotype, genRef, fit), EPSILON);
        // TODO: Test the linear combination
    }

    [Test]
    public void TestReferenceFitness([Values] bool sexXX, [Values(0,1)] int refId, [Values(-1, 0, 1)] int myInt)
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
        var contigs = genRef.GetGenotype(false).Select(region => new Contig(region)).ToList();
        foreach (string chrNo in genRef.AllChrs)
        {
            int chrToCont = chrNo != "chrY" ? genRef.AllChrs.FindIndex(c => c == chrNo) : 45;
            int contigCount = contigs[chrToCont].GetPresentGenes(selectList).Count;
            int chrCount = selectList[chrNo].Count;
            Assert.AreEqual(chrCount, contigCount);
        }
    }
    
    [Test]
    public void TestTsgOgSum([Values] bool sexXX, [Values] bool useTSG, [Values(0,1)] int refId)
    {
        var genRef = _refs[refId];
        var selectList = genRef.GeneLists[useTSG ? GeneListType.TumorSuppressor : GeneListType.Oncogene];
        double sumHap1 = selectList.Where(pair => pair.Key != "chrY").Sum(pair => pair.Value.Sum(g => g.DeltaFitness));
        string missingChr = sexXX ? "chrY" : "chrX";
        double sumHap2 = selectList.Where(pair => pair.Key != missingChr).Sum(pair => pair.Value.Sum(g => g.DeltaFitness));
        double total = sumHap1 + sumHap2;
        var karyotype = new Karyotype(genRef, sexXX);
        karyotype.ApplyWGD();
        var cnList = Fitness.CalcCNs(selectList, karyotype);
        double sum = Fitness.TsgOgTerm(genRef, cnList, sexXX);
        Assert.AreEqual(total, sum, EPSILON);
    }
}
