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
    
    private static Gene MakeGene(ChrNo chrNo, double deltaFitness)
        => new($"G{chrNo}", new Region(0, 50, new ChrID(chrNo, false)), deltaFitness);
    
    private Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> dict = 
        new Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>>();

    [SetUp]
    public void Setup()
    {   
        dict = Enum.GetValues(typeof(GeneListType)).Cast<GeneListType>().ToDictionary(g => g, _ => 
            Enum.GetValues(typeof(ChrNo)).Cast<ChrNo>().ToDictionary(t => t, _ => new List<Gene>()));
        dict[GeneListType.TumorSuppressor][ChrNo.chr2].Add(MakeGene(ChrNo.chr2, 0.01));
        dict[GeneListType.Oncogene][ChrNo.chr4].Add(MakeGene(ChrNo.chr4, 0.01));
        dict[GeneListType.Essentiality][ChrNo.chr3].Add(MakeGene(ChrNo.chr3, 0.01));
        dict[GeneListType.TumorSuppressor][ChrNo.chr1].Add(MakeGene(ChrNo.chr1, 0.01));
        dict[GeneListType.Oncogene][ChrNo.chr1].Add(MakeGene(ChrNo.chr1, 0.02));
        dict[GeneListType.Essentiality][ChrNo.chr1].Add(MakeGene(ChrNo.chr1, 0.03));
        var fparam = new FitnessParams(0.1, 0.1, 0.1);
        Fitness.SetStartingParams(dict, fparam);

        //Fitness.SetStartingParams()
    }

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
        HGRef.Assembly = GenomeAssembly.hg19;
        var karyotype = new Karyotype(true);
        Assert.AreEqual(Fitness.CalcCNs(GeneListType.Essentiality, karyotype).FirstOrDefault(), (dict[GeneListType.Essentiality][ChrNo.chr1].FirstOrDefault(), 2));
        Assert.AreEqual(Fitness.CalcCNs(GeneListType.TumorSuppressor, karyotype).FirstOrDefault(), (dict[GeneListType.TumorSuppressor][ChrNo.chr1].FirstOrDefault(), 2));
        Assert.AreEqual(Fitness.CalcCNs(GeneListType.Oncogene, karyotype).FirstOrDefault(), (dict[GeneListType.Oncogene][ChrNo.chr1].FirstOrDefault(), 2));
        // Delete chr 1 (contig 0)
        karyotype.ApplyContigDeletion(0);
        Assert.AreEqual(Fitness.CalcCNs(GeneListType.Essentiality, karyotype).FirstOrDefault(), (dict[GeneListType.Essentiality][ChrNo.chr1].FirstOrDefault(), 1));
        Assert.AreEqual(Fitness.CalcCNs(GeneListType.TumorSuppressor, karyotype).FirstOrDefault(), (dict[GeneListType.TumorSuppressor][ChrNo.chr1].FirstOrDefault(), 1));
        Assert.AreEqual(Fitness.CalcCNs(GeneListType.Oncogene, karyotype).FirstOrDefault(), (dict[GeneListType.Oncogene][ChrNo.chr1].FirstOrDefault(), 1));
        // Delete chr 2 (contig 1), fitness values for genes on chr 1 should remain the same.
        karyotype.ApplyContigDeletion(1);
        Assert.AreEqual(Fitness.CalcCNs(GeneListType.Essentiality, karyotype).FirstOrDefault(), (dict[GeneListType.Essentiality][ChrNo.chr1].FirstOrDefault(), 1));
        Assert.AreEqual(Fitness.CalcCNs(GeneListType.TumorSuppressor, karyotype).FirstOrDefault(), (dict[GeneListType.TumorSuppressor][ChrNo.chr1].FirstOrDefault(), 1));
        Assert.AreEqual(Fitness.CalcCNs(GeneListType.Oncogene, karyotype).FirstOrDefault(), (dict[GeneListType.Oncogene][ChrNo.chr1].FirstOrDefault(), 1));
    }

    [Test]
    public void TestCalculate()
    {
        HGRef.Assembly = GenomeAssembly.hg19;
        var karyotype = new Karyotype(true);
        var listGenes = Enum.GetValues(typeof(GeneListType)).Cast<GeneListType>().ToDictionary(
            t => t,
            _ => Enum.GetValues(typeof(ChrNo)).Cast<ChrNo>().ToDictionary(chrNo => chrNo, _ => new List<Gene>()));

        listGenes[GeneListType.Oncogene][ChrNo.chr1].Add(MakeGene(ChrNo.chr1, 0.001));
        Assert.AreEqual(1, Fitness.Calculate(karyotype), EPSILON);

        // TODO: Test the linear combination
    }

    [Test]
    public void TestReferenceFitness([Values] bool sexXX, [Values(GenomeAssembly.hg19, GenomeAssembly.hg38)] GenomeAssembly genomeAssembly, [Values(-1, 0, 1)] int myInt)
    {
        const string dataPath = "./../../../../data";
        HGRef.Assembly = genomeAssembly;
        var geneLists = FileIO.ReadGeneLists(dataPath, genomeAssembly);
        HGRef.Assembly = genomeAssembly;
        var karyotype = new Karyotype(sexXX);
        var tsgCNs = Fitness.CalcCNs(GeneListType.TumorSuppressor, karyotype);
        double tsg = Fitness.TsgOgTerm(tsgCNs, sexXX);
        var ogsCNs = Fitness.CalcCNs(GeneListType.Oncogene, karyotype);
        double og = Fitness.TsgOgTerm(ogsCNs, sexXX);
        Console.WriteLine($"sex: {sexXX}, assembly: {genomeAssembly}, TSG: {tsg}, OG: {og}");
    }
    [Test]
    public void TestGetGeneList()
    {
        var karyotype = new Karyotype(true);
        var geneList = new Dictionary<GeneListType, List<Gene>>();
        geneList.Add(GeneListType.TumorSuppressor, new List<Gene>(dict[GeneListType.TumorSuppressor][ChrNo.chr1]));
        geneList.Add(GeneListType.Oncogene, new List<Gene>(dict[GeneListType.Oncogene][ChrNo.chr1]));
        geneList.Add(GeneListType.Essentiality, new List<Gene>(dict[GeneListType.Essentiality][ChrNo.chr1]));
        Assert.AreEqual(Fitness.GetGeneList(0, 1000, ChrNo.chr1), geneList);
        geneList[GeneListType.TumorSuppressor] = new List<Gene>();
        geneList[GeneListType.Oncogene] = new List<Gene>();
        geneList[GeneListType.Essentiality] = new List<Gene>(dict[GeneListType.Essentiality][ChrNo.chr3]);
        Assert.AreEqual(Fitness.GetGeneList(0, 1000, ChrNo.chr3), geneList); 
    }
}
