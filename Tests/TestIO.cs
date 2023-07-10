// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using System;
using System.Collections.Generic;
using System.IO;
using System.Text.Json;
using System.Linq;
using NUnit.Framework;
using SimChA.IO;
using SimChA.DataTypes;
using SimChA.EventData;
using SimChA.Simulation;
using JsonSerializer = System.Text.Json.JsonSerializer;
using System.Text;

namespace Tests;

[TestFixture]
public class TestIO
{
    public static string DATA_PATH = "./../../../../data/";
    public static string HG_19_PATH => DATA_PATH + "./hg19";
    public static string HG_38_PATH => DATA_PATH + "./hg38";
    
    private GenRef _genRef;
    
    [SetUp]
    public void Setup()
    {
        _genRef = FileIO.GetGenRef(HG_19_PATH);
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
        res = Parsers.ParseSimParams(@"{""EventCount"": 10, ""EventDist"": ""Normal""}");
        Assert.AreEqual(10, res.EventCount);
        Assert.AreEqual(Distribution.Normal, res.EventDist);
        res = Parsers.ParseSimParams(@"{""Signatures"": {""test"" : { ""Prob"": 1 }}}");
        Assert.AreEqual("test", res.Signatures!.First().Key);
        Assert.AreEqual(1, res.Signatures!.First().Value.Prob, 0.000001);
        res = Parsers.ParseSimParams(@"{""Signatures"": {""test"" : {""Prob"": 1, ""Events"": [{""Type"": ""WholeGenomeDoubling"", ""Prob"": 0.1}]}}}");
        Assert.AreEqual(CNEventType.WholeGenomeDoubling, res.Signatures!.First().Value.Events.First().Type);
        Assert.AreEqual(0.1, res.Signatures!.First().Value.Events.First().Prob, TestFitness.EPSILON);
    }

    [Test]
    public void TestReadGeneLists()
    {
        var tsgList = _genRef.AllChrs.ToDictionary(t => t, t => new List<Gene>());
        
        string genesTSG = "chr1\t1\t50\tTSG1\t0.001";
        
        var gene1 = new Gene("TSG1", new GenRange(0, 50, "chr1"), 0.001);
        tsgList["chr1"].Add(gene1);
        var listFromString = Parsers.ParseGeneList(new StringReader(genesTSG), _genRef.AllChrs);
        Assert.AreEqual(tsgList, listFromString);
        
        genesTSG += "\nchr2\t100\t5000\tTSG2\t0.01";
        var gene2 = new Gene("TSG2", new GenRange(99, 5000, "chr2"), 0.01);
        tsgList["chr2"].Add(gene2);
        listFromString = Parsers.ParseGeneList(new StringReader(genesTSG), _genRef.AllChrs);
        Assert.AreEqual(tsgList, listFromString);
    }

    [Test]
    public void TestWriteClones()
    {
        // TODO: TestWriteClones does not actually have a test!
        string? projectPath = Path.GetDirectoryName(Path.GetDirectoryName(Path.GetDirectoryName(TestContext.CurrentContext.TestDirectory)));
        var files = new FileIO(projectPath + "/out");
        var kar = new Karyotype( _genRef, false);
        var rnd = new Random(48);
        TestKaryotype.ApplyRandomEvent(rnd, kar, new CNEventPars(CNEventType.Rigma, 1.0, 1_000_000, 10));
        var clone = new CloneIn(1, -1, 0, 1);
    }

    [Test]
    public void TestWriteVCF()
    {
        string? projectPath = Path.GetDirectoryName(Path.GetDirectoryName(Path.GetDirectoryName(TestContext.CurrentContext.TestDirectory)));
        var files = new FileIO(projectPath + "/out");
        const string sequence = 
            @"ACTGACTGACTGACTG";
        const string fasta = @$">chr1
{sequence}";
        byte[] byteArray = Encoding.UTF8.GetBytes(fasta);
        MemoryStream stream = new MemoryStream(byteArray);
        var genContents = Parsers.ParseFasta(new StreamReader(stream)).ToList();
        _genRef.GenContentsDict = new Dictionary<string, StringBuilder> (){{"chr1", genContents[0]}};

        var eventPars = new List<CNEventPars>(){new CNEventPars(CNEventType.SNV, 1.0)};
        
        var clonesIn = new List<CloneIn> (){new CloneIn(0,-1, 0, 1), new CloneIn(1, 0, 1, 1)};

        var sample = new Sample("sample", false, clonesIn, eventPars, null);
        var contigs = new List<Contig> (){new Contig(new Region(0, sequence.Count(), "chr1", true))};
        sample.EventDescs[0] = new List<CNEventDesc> ();
        sample.Kars[0] = new Karyotype(contigs, new List<GenRange>(), false);
        sample.EventDescs[1] = new List<CNEventDesc> ();
        sample.Kars[1] = new Karyotype(sample.Kars[0]);

        long loc = 5;
        int contigID = 0;
        var newNucleotide = Nucleotide.N;
        sample.Kars[1].ApplySNV(contigID, loc, newNucleotide);
        var rnd = new Random(48);
        var eventData = new PointMutationData(rnd, eventPars[0], contigID, sequence.Count());
        var eventDesc = new CNEventDesc(CNEventType.SNV, 1, eventData.ToString());
        sample.EventDescs[1].Append(eventDesc);

        var samples = new List<Sample> (){sample};
        if (samples.Any(s => s.EventDescs.Any()) )
        {
            Console.WriteLine("Hello");
            files.WriteVCF(_genRef, samples);
        }
    }

    [Test]
    public void TestReadClones()
    {
        const string clonesStr = "ID,ParentID,Distance,Alive,Necrotic,Lost,Drivers,Passengers,Fitness,CCF\n" +
                                 "0,0,0,116515,0,584410,1,0,1.2,0.9755\n" +
                                 "1,0,1,149897,0,339690,2,0,1.44,0.8619\n" +
                                 "5,1,1,310270,0,426497,3,0,1.728,0.3025\n" +
                                 "6,1,1,423957,0,583948,3,0,1.728,0.4133";
        var clones = Parsers.ParseClones(new StringReader(clonesStr), true, ",");
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
        Assert.AreEqual(2, profiles["1"].FindRegionsOfChr("chr1").Count()); // 2 existing
        Assert.AreEqual(4, profiles["1"].FindRegionsOfChr("chr2").Count()); // 4 missing (split by null regions)
        Assert.AreEqual(9, profiles["1"].FindRegionsOfChr("chr3").Count()); // 5 existing + 4 missing
        Assert.AreEqual(2, profiles["1"].FindRegionsOfChr("chr4").Count()); // 2 missing
        Assert.AreEqual(false, profiles["2"].SexXX);
    }
    
    [Test]
    public void TestParseChromFile()
    {
        Assert.AreEqual("hg19", _genRef.Name);
        Assert.AreEqual(22, _genRef.AutosomeCount);
        Assert.AreEqual(46, _genRef.ChrCount);
    }

    [Test]
    public void TestParseFasta()
    {
        const string sequence = 
            @"ACTGACTGACTGACTG";
        const string fasta = @$">chr1
{sequence}";
        byte[] byteArray = Encoding.UTF8.GetBytes(fasta);
        MemoryStream stream = new MemoryStream(byteArray);
        var genContents = Parsers.ParseFasta(new StreamReader(stream)).ToList();
        Assert.AreEqual(1, genContents.Count);
        Assert.AreEqual(sequence, genContents[0].ToString());
    }

    [Test]
    public void TestWriteFasta()
    {
        string? projectPath = Path.GetDirectoryName(Path.GetDirectoryName(Path.GetDirectoryName(TestContext.CurrentContext.TestDirectory)));
        var files = new FileIO(projectPath + "/out");
        const string sequence1 = 
            @"AAAAACTGGGGG";
        const string sequence2 =
            @"TTTTTTTTTTTTT";
        const string fasta = @$">chr1
{sequence1}
>chr2
{sequence2}";

        byte[] byteArray = Encoding.UTF8.GetBytes(fasta);
        MemoryStream stream = new MemoryStream(byteArray);
        var genContents = Parsers.ParseFasta(new StreamReader(stream)).ToList();
        _genRef.GenContentsDict = new Dictionary<string, StringBuilder> (){{"chr1", genContents[0]},{"chr2", genContents[1]}};

        var eventPars = new List<CNEventPars>(){new CNEventPars(CNEventType.InternalInversion, 1, 10)};
        
        var clonesIn = new List<CloneIn> (){new CloneIn(0,-1, 1, 1)};

        var sample = new Sample("sample_1", false, clonesIn, eventPars, null);
        var contigs = new List<Contig> (){new Contig(new Region(0, sequence1.Count(), "chr1", true)), new Contig(new Region(0, sequence2.Count(), "chr2", true))};
        sample.EventDescs[0] = new List<CNEventDesc> ();
        sample.Kars[0] = new Karyotype( contigs, new List<GenRange>(), false);
        // Apply an internal inversion
        var rnd = new Random(0);
        var contigID = 0;
        var eventData = new InternalEventData(rnd, eventPars[0], contigID, sequence1.Count());
        var eventDesc = new CNEventDesc(CNEventType.InternalInversion, 1, eventData.ToString());
        sample.EventDescs[0].Append(eventDesc);
        sample.Kars[0].ApplyInternalInversion(contigID, 4, 8);

        contigID = 1;
        var snvEventPars = new CNEventPars(CNEventType.SNV, 1);
        var snvEventData = new PointMutationData(rnd, snvEventPars, contigID, sequence2.Count());
        var snvEventDesc = new CNEventDesc(CNEventType.SNV, 2, snvEventData.ToString());
        sample.EventDescs[0].Append(snvEventDesc);
        sample.Kars[0].ApplySNV(contigID, 10, Nucleotide.A);

        var samples = new List<Sample> (){sample};
        if (samples.Any(s => s.EventDescs.Any()) )
        {
            files.WriteFasta(_genRef, samples);
        }
        // The expected output is
        /*
        >ctg0
        AAAAGTCAGGGG
        >ctg1
        TTTTTTTTTTATT
        */
    }
}