// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using System;
using System.Collections.Generic;
using System.IO;
using System.Text.Json;
using System.Linq;
using SimChA.IO;
using SimChA.DataTypes;
using SimChA.EventData;
using SimChA.Simulation;
using NUnit.Framework;
using JsonSerializer = System.Text.Json.JsonSerializer;

namespace Tests;

[TestFixture]
public class TestIO
{
    private GenRef _genRef;
    
    [SetUp]
    public void Setup()
    {
        _genRef = Parsers.ParseChromosomes("test", TestData.TEST_CHROMOSOMES);
    }

    
    [Test]
    public void TestContig()
    {
        var clone = new CloneIn(-1, 0, 1, 1);
        Assert.DoesNotThrow(() => clone.ToString());
    }
    
    [Test]
    public void TestConfigSerialization()
    {
        var fit = new FitnessParams(0.001f, 0.01f, 0.000_1f);
        var simParams = new SimParams(0, SexEnum.Both, 1, Distribution.Uniform, fit, null, null);
        var options = new JsonSerializerOptions { WriteIndented = true };
        string serialized = JsonSerializer.Serialize(simParams, options);
        Console.WriteLine(serialized);
        var deserialized = JsonSerializer.Deserialize<SimParams>(serialized);
        Assert.NotNull(deserialized);
        Assert.AreEqual(simParams, deserialized);
    }

    [Test]
    public void TestConfigRead()
    {
        var res = Parsers.ParseSimParams(@"{}");
        Assert.AreEqual(0, res.Seed);
        res = Parsers.ParseSimParams(@"{""EventCount"": 10, ""Distribution"": ""Normal""}");
        Assert.AreEqual(10, res.EventCount);
        Assert.AreEqual(Distribution.Normal, res.Distribution);
        res = Parsers.ParseSimParams(@"{""Signatures"": [{""Name"": ""test"", ""Prob"": 1}]}");
        Assert.AreEqual(1, res.Signatures!.First().Prob, 0.000001);
        res = Parsers.ParseSimParams(@"{""Signatures"": [{""Name"": ""test"", ""Prob"": 1, ""Events"": [{""Type"": ""WholeGenomeDoubling"", ""Prob"": 0.1}]}]}");
        Assert.AreEqual(CNEventType.WholeGenomeDoubling, res.Signatures!.First().Events.First().Type);
        res = Parsers.ParseSimParams(@"{""Signatures"": [{""Name"": ""test"", ""Prob"": 1, ""Events"": [{""Type"": ""InternalInversion"", ""Prob"": 0.1, ""Pars"": {""Mean"": 0.1}}]}]}");
        Assert.AreEqual(0.1, res.Signatures!.First().Events.First().Pars!["Mean"], 0.000001);
    }

    [Test]
    public void TestReadGeneLists()
    {
        var tsgList = Enum.GetValues<ChrNo>().ToDictionary(t => t, t => new List<Gene>());
        
        string genesTSG = "chr1\t1\t50\tTSG1\t0.001";
        
        var gene1 = new Gene("TSG1", new GenRange(0, 50, ChrNo.chr1), 0.001);
        tsgList[ChrNo.chr1].Add(gene1);
        var listFromString = Parsers.ParseGeneList(new StringReader(genesTSG));
        Assert.AreEqual(tsgList, listFromString);
        
        genesTSG += "\nchr2\t100\t5000\tTSG2\t0.01";
        var gene2 = new Gene("TSG2", new GenRange(99, 5000, ChrNo.chr2), 0.01);
        tsgList[ChrNo.chr2].Add(gene2);
        listFromString = Parsers.ParseGeneList(new StringReader(genesTSG));
        Assert.AreEqual(tsgList, listFromString);
    }

    [Test]
    public void TestWriteClones()
    {
        string? projectPath = Path.GetDirectoryName(Path.GetDirectoryName(Path.GetDirectoryName(TestContext.CurrentContext.TestDirectory)));
        var files = new FileIO(projectPath + "/out");
        var kar = new Karyotype( _genRef, false);
        var pars = new Dictionary<string, double> { ["Size"] = 1_000_000, ["Frag"] = 10 };
        var rnd = new Random(48);
        TestKaryotype.ApplyRandomEvent(rnd, kar, new CNEventPars(CNEventType.Rigma, 1.0, pars));
        var clone = new CloneIn(1, -1, 0, 1);
    }

    [Test]
    public void TestReadClones()
    {
        const string clonesStr = "ID,ParentID,Distance,Alive,Necrotic,Lost,Drivers,Passengers,Fitness,CCF\n" +
                                 "0,0,0,116515,0,584410,1,0,1.2,0.9755\n" +
                                 "1,0,1,149897,0,339690,2,0,1.44,0.8619\n" +
                                 "5,1,1,310270,0,426497,3,0,1.728,0.3025\n" +
                                 "6,1,1,423957,0,583948,3,0,1.728,0.4133";
        var clones = Parsers.ParseClones(new StringReader(clonesStr), true);
        Assert.AreEqual(4, clones.Count);
        Assert.AreEqual(0, clones[0].CloneId);
        Assert.AreEqual(0, clones[1].ParentId);
        Assert.AreEqual(1, clones[2].Distance);
        Assert.AreEqual(1.728, clones[3].FitnessTarget, double.Epsilon * 10);
        
    }

    [Test]
    public void TestParseCNAProfiles()
    {
        const string Profiles = @"sample_id	chrom	start	end	cn_a	cn_b
1	chr1	1	249250621	1	1
1	chr2	13133	2429856	0	0
1	chr3	62226	171636043	2	3
2	chrX	2	6	1	0
2	chrY	3	4	0	1";
        var profiles = Parsers.ParseCNAProfile(_genRef, new StringReader(Profiles));
        Assert.AreEqual(2, profiles.Count);
        Assert.AreEqual(2, profiles["1"].CountContigs());
        Assert.AreEqual(2, profiles["1"].FindRegionsOfChr(ChrNo.chr1).Count()); // 2 existing
        Assert.AreEqual(4, profiles["1"].FindRegionsOfChr(ChrNo.chr2).Count()); // 4 missing (split by null regions)
        Assert.AreEqual(9, profiles["1"].FindRegionsOfChr(ChrNo.chr3).Count()); // 5 existing + 4 missing
        Assert.AreEqual(2, profiles["1"].FindRegionsOfChr(ChrNo.chr4).Count()); // 2 missing
        Assert.AreEqual(false, profiles["2"].SexXX);
    }

    [Test]
    public void TestParseChromText()
    {
        const string testRef = @"chr1	249250621
chr2	243199373
chr3	198022430	Both
chrX	155270560	0
chrY	59373566";
        var genRef = Parsers.ParseChromosomes("testRef", testRef);
        Assert.AreEqual(3, genRef.AutosomeCount);
        Assert.AreEqual(8, genRef.ChrCount);
        Assert.AreEqual(249250621, genRef.GetChromLen(ChrNo.chr1));
        Assert.AreEqual(198022430, genRef.GetChromLen(ChrNo.chr3));
        Assert.AreEqual(155270560, genRef.GetChromLen(ChrNo.chrX));
        Assert.AreEqual(59373566, genRef.GetChromLen(ChrNo.chrY));
    }

    [Test]
    public void TestParseChromFile()
    {
        const string dataPath = "./../../../../data/hg19";
        var genRef = FileIO.ReadChromosomes(dataPath);
        Assert.AreEqual("hg19", genRef.Name);
        Assert.AreEqual(22, genRef.AutosomeCount);
        Assert.AreEqual(46, genRef.ChrCount);
    }
}