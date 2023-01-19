// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Simulation;

namespace Tests;

[TestFixture]
public class TestFitness
{
    private static Gene MakeGene(ChrNo chrNo, double deltaFitness)
        => new("G" + chrNo, new Region(0, 50, new ChrID(chrNo, false)), deltaFitness);

    [Test]
    public void TestEssTerm()
    {
        Assert.AreEqual(0, Fitness.EssTerm(new List<(Gene, int)>()));
        
        var testNoEffect = new List<(Gene, int)> { (MakeGene(ChrNo.chr1, 0), 0) };
        Assert.AreEqual(0, Fitness.EssTerm(testNoEffect));
        
        var testMissing = new List<(Gene, int)> { (MakeGene(ChrNo.chr1, 0.1), 0) };
        Assert.AreEqual(0.1, Fitness.EssTerm(testMissing));

        var testHaplosufficient = new List<(Gene, int)> { (MakeGene(ChrNo.chr1, 0.1), 1) };
        Assert.AreEqual(0, Fitness.EssTerm(testHaplosufficient));
        
        var testList = new List<(Gene, int)> { (MakeGene(ChrNo.chr1, 0.1), 0), (MakeGene(ChrNo.chr2, 0.2), 0) };
        Assert.AreEqual(0.1 + 0.2, Fitness.EssTerm(testList));
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
        Assert.AreEqual(0, Fitness.StressTerm(46));
        Assert.AreEqual(0, Fitness.StressTerm(45));
        Assert.AreEqual(1, Fitness.StressTerm(47));
        Assert.AreEqual(4, Fitness.StressTerm(48));
    }
    
    [Test]
    public void TestCNCalulation()
    {
        // Seed 14 to get chr1 delete
        var rnd = new Random(14);
        var karyotype = new Karyotype(true);
        var dict = Enum.GetValues(typeof(ChrNo)).Cast<ChrNo>().ToDictionary(t => t, t => new List<Gene>());
        dict[ChrNo.chr1].Add(MakeGene(ChrNo.chr1, 0.01));
        Assert.AreEqual(Fitness.CalcCNs(dict, karyotype).FirstOrDefault(), (dict[ChrNo.chr1].FirstOrDefault(), 2));
        karyotype.ApplyAberration(rnd, AberrationEnum.ChromDeletion, new SimChA.IO.BaseAbbP(1));
        Assert.AreEqual(Fitness.CalcCNs(dict, karyotype).FirstOrDefault(), (dict[ChrNo.chr1].FirstOrDefault(), 1));
        karyotype.ApplyAberration(rnd, AberrationEnum.ChromDeletion, new SimChA.IO.BaseAbbP(1));
        Assert.AreEqual(Fitness.CalcCNs(dict, karyotype).FirstOrDefault(), (dict[ChrNo.chr1].FirstOrDefault(), 1));
    }

    [Test]
    public void TestCalculate()
    {
        var karyotype = new Karyotype(true);
        var simParams = SimParams.CreateSimParams(14, true, 0.001f, 0.01f, 0.000_1f, AberrationsInfo.DefaultAberrations());
        //var listGenes = new List<Dictionary<ChrNo, List<Gene>>>();
        var listGenes = Enum.GetValues(typeof(GeneListType)).Cast<GeneListType>().ToDictionary(
            t => t, 
            t => Enum.GetValues(typeof(ChrNo)).Cast<ChrNo>().ToDictionary(t => t, t => new List<Gene>()));
        
        
        listGenes[GeneListType.Oncogene][ChrNo.chr1].Add(MakeGene(ChrNo.chr1, 0.001));
        Assert.AreEqual(1, Fitness.Calculate(karyotype, listGenes, simParams));

        //For OGs with one chromosome lost
        karyotype.ApplyAberration(new Random(14), AberrationEnum.ChromDeletion, new BaseAbbP(1));
        Assert.AreEqual(1+(1-2)*(simParams.TsgOgFraction*0.001), Fitness.Calculate(karyotype, listGenes, simParams), 0.0000000001);
        
        //For TSGs with one chromosome lost
        listGenes[GeneListType.Oncogene][ChrNo.chr1].Clear();
        listGenes[GeneListType.TumorSuppressor][ChrNo.chr1].Add(MakeGene(ChrNo.chr1, 0.0001));
        Assert.AreEqual(1+(1-2)*(-simParams.TsgOgFraction*0.0001), Fitness.Calculate(karyotype, listGenes, simParams),  0.0000000001);
        
        // For essential genes with one chromosome lost
        listGenes[GeneListType.TumorSuppressor][ChrNo.chr1].Clear();
        listGenes[GeneListType.Essentiality][ChrNo.chr1].Add(MakeGene(ChrNo.chr1, 0.01));
        Assert.AreEqual(1, Fitness.Calculate(karyotype, listGenes, simParams));
        
        //Seed to lose second chromosome 1
        karyotype.ApplyAberration(new Random(77), AberrationEnum.ChromDeletion, new BaseAbbP(1));
        
        //For OGs with chromosome 1 lost twice
        listGenes[GeneListType.Essentiality][ChrNo.chr1].Clear();
        listGenes[GeneListType.Oncogene][ChrNo.chr1].Add(MakeGene(ChrNo.chr1, 0.001));
        Assert.AreEqual(1+(0-2)*(simParams.TsgOgFraction*0.001), Fitness.Calculate(karyotype, listGenes, simParams), 0.0000000001);
        
        //For TSGs with chromosome 1 lost twice
        listGenes[GeneListType.Oncogene][ChrNo.chr1].Clear();
        listGenes[GeneListType.TumorSuppressor][ChrNo.chr1].Add(MakeGene(ChrNo.chr1, 0.0001));
        Assert.AreEqual(1+(0-2)*(-simParams.TsgOgFraction*0.0001), Fitness.Calculate(karyotype, listGenes, simParams), 0.0000000001);
        
        //For essential genes with chromosome 1 lost twice
        listGenes[GeneListType.TumorSuppressor][ChrNo.chr1].Clear();
        listGenes[GeneListType.Essentiality][ChrNo.chr1].Add(MakeGene(ChrNo.chr1, 0.01));
        Assert.AreEqual(1+(-1)*(simParams.EssentialFraction*0.01), Fitness.Calculate(karyotype, listGenes, simParams), 0.0000000001);
    }
}