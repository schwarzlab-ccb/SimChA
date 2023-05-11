// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
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
    
    private static Gene MakeGene(ChrNo chrNo, double deltaFitness)
        => new($"G{chrNo}", new Region(0, 50, new ChrID(chrNo, false)), deltaFitness);

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
        Assert.AreEqual(0, Fitness.TsgOgTerm(new List<(Gene, int)>()));
        
        var testNoEffect = new List<(Gene, int)> {(MakeGene(ChrNo.chr1, 0), 0)};
        Assert.AreEqual(0, Fitness.TsgOgTerm(testNoEffect));

        var testMissing = new List<(Gene, int)> { (MakeGene(ChrNo.chr1, 0.1), 1)};
        Assert.AreEqual(-0.1, Fitness.TsgOgTerm(testMissing));

        var testMissingTwice = new List<(Gene, int)> {(MakeGene(ChrNo.chr1, 0.1), 0)};
        Assert.AreEqual(-0.2, Fitness.TsgOgTerm(testMissingTwice));

        var testList = new List<(Gene, int)> {
            (MakeGene(ChrNo.chr1, 0.1), 1), 
            (MakeGene(ChrNo.chr1, 0.2), 0), 
            (MakeGene(ChrNo.chr1, 0.3), 2)
        };
        Assert.AreEqual(-0.2 - 0.2 - 0.1, Fitness.TsgOgTerm(testList));
    }
    
    [Test]
    public void TestStressTerm()
    {
        var xxKaryotype = new Karyotype(true);
        var xyKaryotype = new Karyotype(false);
        Assert.AreEqual(0, Fitness.StressTerm(xxKaryotype.GenomeLen(), true), EPSILON);
        Assert.AreEqual(0, Fitness.StressTerm(xyKaryotype.GenomeLen(), false), EPSILON);
        Assert.AreNotEqual(0, Fitness.StressTerm(xxKaryotype.GenomeLen(), false));
        xxKaryotype.ApplyWGD(); // Double all
        Assert.AreEqual(-1, Fitness.StressTerm(xxKaryotype.GenomeLen(), true), EPSILON);
        foreach (int i in Enumerable.Range(0, HGRef.CHR_COUNT)) { xyKaryotype.ApplyContigDeletion(i); }
        Assert.AreEqual(1, Fitness.StressTerm(xyKaryotype.GenomeLen(), false), EPSILON);
    }
    
    [Test]
    public void TestCNCalulation()
    {
        // Seed 14 to get chr1 delete
        var rnd = new Random(14);
        var karyotype = new Karyotype(true);
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
    public void TestCalculate()
    {
        var karyotype = new Karyotype(true);
        var fit = new FitnessParams(0.001f, 0.01f, 0.000_1f);
        //var listGenes = new List<Dictionary<ChrNo, List<Gene>>>();
        var listGenes = Enum.GetValues(typeof(GeneListType)).Cast<GeneListType>().ToDictionary(
            t => t,
            _ => Enum.GetValues(typeof(ChrNo)).Cast<ChrNo>().ToDictionary(chrNo => chrNo, _ => new List<Gene>()));

        listGenes[GeneListType.Oncogene][ChrNo.chr1].Add(MakeGene(ChrNo.chr1, 0.001));
        Assert.AreEqual(1, Fitness.Calculate(karyotype, listGenes, fit), EPSILON);

        // TODO: Test the linear combination
    }

    [Test]
    public void TestOnDataMale()
    {
        const bool sexXX = false;
        const string dataPath = "./../../../../data";
        var geneLists = FileIO.ReadGeneLists(dataPath, sexXX, GenomeAssembly.hg19);
        var karyotype = new Karyotype(sexXX);
        var fit = new FitnessParams(1, 1, 1);
        var fitness = Fitness.Calculate(karyotype, geneLists, fit);
        Assert.AreEqual(1.0, fitness, EPSILON);
    }

    [Test]
    public void TestOnDataFemale()
    {
        const bool sexXX = true;
        const string dataPath = "./../../../../data";
        var geneLists = FileIO.ReadGeneLists(dataPath, sexXX, GenomeAssembly.hg19);
        var karyotype = new Karyotype(sexXX);
        var fit = new FitnessParams(1, 1, 1);
        var fitness = Fitness.Calculate(karyotype, geneLists, fit);
        Assert.AreEqual(1.0, fitness, EPSILON);
    }
}