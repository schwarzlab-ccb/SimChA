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
    private const double EPSILON = 0.0000000001;

    private Dictionary<GenomeAssembly, GenRef> _refs;

    [SetUp]
    public void Setup()
    {
        var hg19Path = "./../../../../data/hg19";
        var hg38Path = "./../../../../data/hg38";
        _refs = new Dictionary<GenomeAssembly, GenRef>
        {
            [GenomeAssembly.hg19] = FileIO.ReadChromosomes(hg19Path),
            [GenomeAssembly.hg38] = FileIO.ReadChromosomes(hg38Path)
        };
        _refs[GenomeAssembly.hg19].GeneLists = FileIO.ReadGeneLists(hg19Path);
        _refs[GenomeAssembly.hg38].GeneLists = FileIO.ReadGeneLists(hg38Path);
    }
    
    private static Gene MakeGene(ChrNo chrNo, double deltaFitness)
        => new($"G{chrNo}", new Region(0, 50, chrNo, false), deltaFitness);

    [Test]
    public void TestEssTerm()
    {
        Assert.AreEqual(0, Fitness.EssTerm(new List<(Gene, int)>()));
        
        var testNoEffect = new List<(Gene, int)> { (MakeGene(ChrNo.chr1, 0), 0) };
        Assert.AreEqual(0, Fitness.EssTerm(testNoEffect));
        
        var testMissing = new List<(Gene, int)> { (MakeGene(ChrNo.chr1, 0.1), 0) };
        Assert.AreEqual(-0.1, Fitness.EssTerm(testMissing));

        var testHaplosufficient = new List<(Gene, int)> { (MakeGene(ChrNo.chr1, 0.1), 1) };
        Assert.AreEqual(0, Fitness.EssTerm(testHaplosufficient));
        
        var testList = new List<(Gene, int)> { (MakeGene(ChrNo.chr1, 0.1), 0), (MakeGene(ChrNo.chr2, 0.2), 0) };
        Assert.AreEqual(-0.1 + -0.2, Fitness.EssTerm(testList));
    }

    [Test]
    public void TestTsgOgTerm()
    {
        Assert.AreEqual(0, Fitness.TsgOgTerm(new List<(Gene, int)>(), true));
        
        var testNoEffect = new List<(Gene, int)> {(MakeGene(ChrNo.chr1, 0), 0)};
        Assert.AreEqual(0, Fitness.TsgOgTerm(testNoEffect, true));

        var testMissing = new List<(Gene, int)> { (MakeGene(ChrNo.chr1, 0.1), 1)};
        Assert.AreEqual(-0.1, Fitness.TsgOgTerm(testMissing, true));

        var testMissingTwice = new List<(Gene, int)> {(MakeGene(ChrNo.chr1, 0.1), 0)};
        Assert.AreEqual(-0.2, Fitness.TsgOgTerm(testMissingTwice, true));

        var testList = new List<(Gene, int)> {
            (MakeGene(ChrNo.chr1, 0.1), 1), 
            (MakeGene(ChrNo.chr1, 0.2), 0), 
            (MakeGene(ChrNo.chr1, 0.3), 2)
        };
        Assert.AreEqual(-0.2 - 0.2 - 0.1, Fitness.TsgOgTerm(testList, true));

        testList = new List<(Gene, int)>{
            (MakeGene(ChrNo.chrX, 0.1), 2)
        };
        Assert.AreEqual(0, Fitness.TsgOgTerm(testList, true));
        Assert.AreEqual(0.1, Fitness.TsgOgTerm(testList, false));
        testList = new List<(Gene, int)>{
            (MakeGene(ChrNo.chrX, 0.2), 2),
            (MakeGene(ChrNo.chrY, 0.1), 2)
        };
        Assert.AreEqual(0.2, Fitness.TsgOgTerm(testList, true));
        Assert.AreEqual(0.2 + 0.1, Fitness.TsgOgTerm(testList, false));
    }
    
    [Test]
    public void TestStressTerm([Values(GenomeAssembly.hg19, GenomeAssembly.hg38)] GenomeAssembly genomeAssembly)
    {
        var genRef = _refs[genomeAssembly];
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
    public void TestCNCalulation([Values(GenomeAssembly.hg19, GenomeAssembly.hg38)] GenomeAssembly genomeAssembly)
    {
        // Seed 14 to get chr1 delete
        var genRef = _refs[genomeAssembly];
        var rnd = new Random(14);
        var karyotype = new Karyotype(genRef, true);
        var deletion = new CNEventPars(CNEventType.ChromDeletion, 1);
        var dict = Enum.GetValues(typeof(ChrNo)).Cast<ChrNo>().ToDictionary(t => t, _ => new List<Gene>());
        dict[ChrNo.chr1].Add(MakeGene(ChrNo.chr1, 0.01));
        Assert.AreEqual(Fitness.CalcCNs(dict, karyotype).FirstOrDefault(), (dict[ChrNo.chr1].FirstOrDefault(), 2));
        TestKaryotype.ApplyRandomEvent(rnd, karyotype, deletion);
        Assert.AreEqual(Fitness.CalcCNs(dict, karyotype).FirstOrDefault(), (dict[ChrNo.chr1].FirstOrDefault(), 1));
        TestKaryotype.ApplyRandomEvent(rnd, karyotype, deletion);
        Assert.AreEqual(Fitness.CalcCNs(dict, karyotype).FirstOrDefault(), (dict[ChrNo.chr1].FirstOrDefault(), 1));
    }

    [Test]
    public void TestCalculate([Values(GenomeAssembly.hg19, GenomeAssembly.hg38)] GenomeAssembly genomeAssembly)
    {
        var genRef = _refs[genomeAssembly];
        var karyotype = new Karyotype(genRef, true);
        var fit = new FitnessParams(0.001f, 0.01f, 0.000_1f);
        genRef.GeneLists = Enum.GetValues(typeof(GeneListType)).Cast<GeneListType>().ToDictionary(
            t => t,
            _ => Enum.GetValues(typeof(ChrNo)).Cast<ChrNo>().ToDictionary(chrNo => chrNo, _ => new List<Gene>()));
        genRef.GeneLists[GeneListType.Oncogene][ChrNo.chr1].Add(MakeGene(ChrNo.chr1, 0.001));
        Assert.AreEqual(1, Fitness.Calculate(karyotype, genRef, fit), EPSILON);

        // TODO: Test the linear combination
    }

    [Test]
    public void TestReferenceFitness([Values] bool sexXX, [Values(GenomeAssembly.hg19, GenomeAssembly.hg38)] GenomeAssembly genomeAssembly, [Values(-1, 0, 1)] int myInt)
    {
        var genRef = _refs[genomeAssembly];
        var karyotype = new Karyotype(genRef, sexXX);
        var tsgCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.TumorSuppressor], karyotype);
        double tsg = Fitness.TsgOgTerm(tsgCNs, sexXX);
        Assert.AreEqual(0, tsg, EPSILON);
        var ogsCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.Oncogene], karyotype);
        double og = Fitness.TsgOgTerm(ogsCNs, sexXX);;
        Assert.AreEqual(0, og, EPSILON);
    }


    [Test]
    public void TestGetPresentGenes([Values] bool useTSG, [Values(GenomeAssembly.hg19, GenomeAssembly.hg38)] GenomeAssembly genomeAssembly)
    {
        var genRef = _refs[genomeAssembly];
        var selectList = genRef.GeneLists[useTSG ? GeneListType.TumorSuppressor : GeneListType.Oncogene];
        var contigs = genRef.GetGenotype(false).Select(region => new Contig(region)).ToList();
        foreach (ChrNo chrNo in Enum.GetValues(typeof(ChrNo)))
        {
            int chrToCont = chrNo != ChrNo.chrY ? (int)chrNo : 45;
            int contigCount = contigs[chrToCont].GetPresentGenes(selectList).Count;
            int chrCount = selectList[chrNo].Count;
            Assert.AreEqual(chrCount, contigCount);
        }
    }
    
    [Test]
    public void TestTsgOgSum([Values] bool sexXX, [Values] bool useTSG, [Values(GenomeAssembly.hg19, GenomeAssembly.hg38)] GenomeAssembly genomeAssembly)
    {
        var genRef = _refs[genomeAssembly];
        var selectList = genRef.GeneLists[useTSG ? GeneListType.TumorSuppressor : GeneListType.Oncogene];
        double sumHap1 = selectList.Where(pair => pair.Key != ChrNo.chrY).Sum(pair => pair.Value.Sum(g => g.DeltaFitness));
        var missingChr = sexXX ? ChrNo.chrY : ChrNo.chrX;
        double sumHap2 = selectList.Where(pair => pair.Key != missingChr).Sum(pair => pair.Value.Sum(g => g.DeltaFitness));
        double total = sumHap1 + sumHap2;
        var karyotype = new Karyotype(genRef, sexXX);
        karyotype.ApplyWGD();
        var cnList = Fitness.CalcCNs(selectList, karyotype);
        double sum = Fitness.TsgOgTerm(cnList, sexXX);
        Assert.AreEqual(total, sum, EPSILON);
    }
}