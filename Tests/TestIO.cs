// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using NUnit.Framework;
using SimChA.DataTypes;
using SimChA.EventData;
using SimChA.IO;
using SimChA.Simulation;
using JsonSerializer = System.Text.Json.JsonSerializer;

namespace Tests;

[TestFixture]
public class TestIO
{
    public const string DATA_PATH = "./../../../../data/";
    public const string HG_19_PATH = DATA_PATH + "./hg19";
    public const string HG_38_PATH = DATA_PATH + "./hg38";
    private GenRef _genRef;

    [SetUp]
    public void Setup()
    {
        _genRef = FileIO.GetGenRef(HG_19_PATH);
    }

    [Test]
    public void TestRef()
    {
        Console.WriteLine(_genRef);
    }

    [Test]
    public void TestContig()
    {
        var clone = new CloneIn("-1", "0", 1, 1);
        Assert.DoesNotThrow(() => clone.ToString());
    }

    [Test]
    public void TestConfigSerialization()
    {
        var fit = new FitnessParams(0.001f, 0.01f, 0.000_1f, 1f);
        var autosomesOnly = false;
        var simParams = new SimParams(0, SexEnum.None, autosomesOnly, 1, Distribution.Uniform, fit);
        var options = new JsonSerializerOptions { WriteIndented = true };
        string serialized = JsonSerializer.Serialize(simParams, options);
        var deserialized = JsonSerializer.Deserialize<SimParams>(serialized);
        Assert.NotNull(deserialized);
        // Assure that the deserialized object is the same as the original, including the nested objects
        Assert.AreEqual(simParams.Seed, deserialized!.Seed);
        Assert.AreEqual(simParams.Fitness.TotalStrength, deserialized.Fitness.TotalStrength);
        Assert.AreEqual(simParams.Signatures, deserialized.Signatures);
    }

    [Test]
    public void TestConfigRead()
    {
        var res = Parsers.ParseSimParams(@"{}");
        Assert.AreEqual(0, res.Seed);
        res = Parsers.ParseSimParams(@"{""EventCountMean"": 10, ""EventDist"": ""Normal""}");
        Assert.AreEqual(10, res.EventCountMean);
        Assert.AreEqual(Distribution.Normal, res.EventDist);
        res = Parsers.ParseSimParams(@"{""Signatures"": {""test"" : { ""Prob"": 1 }}}");
        Assert.AreEqual("test", res.Signatures!.First().Key);
        Assert.AreEqual(1, res.Signatures!.First().Value.Prob, 0.000001);
        res = Parsers.ParseSimParams(
            @"{""Signatures"": {""test"" : {""Prob"": 1, ""Events"": [{""Type"": ""WholeGenomeDoubling"", ""Prob"": 0.1}]}}}");
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
        string? projectPath =
            Path.GetDirectoryName(Path.GetDirectoryName(Path.GetDirectoryName(TestContext.CurrentContext.TestDirectory)));
        var files = new FileIO(projectPath + "/out");
        var kar = new Karyotype(_genRef, SexEnum.Male);
        var rnd = new Random(48);
        TestKaryotype.ApplyRandomEvent(rnd, kar, new CNEventPars(CNEventType.Rigma, 1.0, 1_000_000, 10));
        var clone = new CloneIn("1", "-1", 0, 1);
    }

    [Test]
    public void TestWriteVCF()
    {
        string? projectPath =
            Path.GetDirectoryName(Path.GetDirectoryName(Path.GetDirectoryName(TestContext.CurrentContext.TestDirectory)));
        var files = new FileIO(projectPath + "/out");
        const string sequence = "ACTGACTGACTGACTG";
        const string fasta = $">chr1\n{sequence}";
        byte[] byteArray = Encoding.UTF8.GetBytes(fasta);
        var stream = new MemoryStream(byteArray);
        var genContents = Parsers.ParseFasta(new StreamReader(stream)).ToList();
        _genRef.GenContentsDict = new Dictionary<string, StringBuilder> { { "chr1", genContents[0] } };

        var eventPars = new List<CNEventPars>
            { new(CNEventType.SNV, 1.0) };

        var clonesIn = new List<CloneIn>
            { new("0", "-1", 0, 1), new("1", "0", 1, 1) };

        var sample = new Sample("sample", SexEnum.Male, clonesIn, eventPars, null);
        var contigs = new List<Contig> { new(new Region(0, sequence.Length, "chr1", true)) };
        sample.EventDescs["0"] = new List<CNEventDesc>();
        sample.Kars["0"] = new Karyotype(contigs, new List<GenRange>(), _genRef.Centromeres, SexEnum.Male);
        sample.EventDescs["1"] = new List<CNEventDesc>();
        sample.Kars["1"] = new Karyotype(sample.Kars["0"]);

        long loc = 5;
        int contigID = 0;
        var newNucleotide = Nucleotide.N;
        sample.Kars["1"].ApplySNV(contigID, loc, newNucleotide);
        var rnd = new Random(48);
        var eventData = new PointMutationData(rnd, eventPars[0], contigID, sequence.Length);
        var eventDesc = new CNEventDesc(CNEventType.SNV, 1, eventData.ToString());
        sample.EventDescs["1"].Add(eventDesc);
        var samples = new List<Sample> { sample };
        if (samples.Any(s => s.EventDescs.Any()))
        {
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
        var clones = Parsers.ParseClonesWithEvents(new StringReader(clonesStr), true, ",");
        Assert.AreEqual(4, clones.Count);
        Assert.AreEqual("0", clones[0].CloneId);
        Assert.AreEqual("0", clones[1].ParentId);
        Assert.AreEqual(1, clones[2].Distance);
        Assert.AreEqual(1.728, clones[3].FitnessTarget, double.Epsilon * 10);

    }

    [Test]
    public void TestParseCNAProfiles()
    {
        const string dummyFile = "sample_id\tchrom\tstart\tend\tcn_a\tcn_b\n" +
                                "1\tchr1\t1\t249250621\t1\t1\n" +
                                "1\tchr2\t13133\t2429856\t0\t0\n" +
                                "1\tchr3\t62226\t171636043\t2\t3\n" +
                                "2\tchrX\t2\t6\t1\t0\n" +
                                "2\tchrY\t3\t4\t0\t1";
        var profiles = Parsers.ParseCNAProfile(_genRef, new StringReader(dummyFile), false);
        Assert.AreEqual(2, profiles.Count);
        Assert.AreEqual(2, profiles["1"].CountContigs());
        Assert.AreEqual(2, profiles["1"].FindRegionsOfChr("chr1").Count()); // 2 existing
        Assert.AreEqual(4, profiles["1"].FindRegionsOfChr("chr2").Count()); // 4 missing (split by null regions)
        Assert.AreEqual(9, profiles["1"].FindRegionsOfChr("chr3").Count()); // 5 existing + 4 missing
        Assert.AreEqual(2, profiles["1"].FindRegionsOfChr("chr4").Count()); // 2 missing
        Assert.AreEqual(SexEnum.Male, profiles["2"].Sex);
    }

    [Test]
    public void TestParseChromFile()
    {
        Assert.AreEqual("hg19", _genRef.Name);
        Assert.AreEqual(22, _genRef.AutosomesCount);
        Assert.AreEqual(46, _genRef.ChrCount(SexEnum.Female, true));
    }

    [Test]
    public void TestParseFasta()
    {
        const string sequence =
            @"ACTGACTGACTGACTG";
        const string fasta = @$">chr1
{sequence}";
        byte[] byteArray = Encoding.UTF8.GetBytes(fasta);
        var stream = new MemoryStream(byteArray);
        var genContents = Parsers.ParseFasta(new StreamReader(stream)).ToList();
        Assert.AreEqual(1, genContents.Count);
        Assert.AreEqual(sequence, genContents[0].ToString());
    }

    [Test]
    public void TestWriteFasta()
    {
        string? projectPath =
            Path.GetDirectoryName(Path.GetDirectoryName(Path.GetDirectoryName(TestContext.CurrentContext.TestDirectory)));
        var files = new FileIO(projectPath + "/out");
        const string sequence1 = "AAAAACTGGGGG";
        const string sequence2 = "TTTTTTTTTTTTT";
        const string fasta = @$">chr1
{sequence1}
>chr2
{sequence2}";

        byte[] byteArray = Encoding.UTF8.GetBytes(fasta);
        var stream = new MemoryStream(byteArray);
        var genContents = Parsers.ParseFasta(new StreamReader(stream)).ToList();
        _genRef.GenContentsDict = new Dictionary<string, StringBuilder>
            { { "chr1", genContents[0] }, { "chr2", genContents[1] } };

        var eventPars = new List<CNEventPars>
            { new(CNEventType.InternalInversion, 1, 10) };

        var clonesIn = new List<CloneIn>
            { new("0", "-1", 1, 1) };

        var sample = new Sample("sample_1", SexEnum.Male, clonesIn, eventPars, null);
        var contigs = new List<Contig>
            { new(new Region(0, sequence1.Count(), "chr1", true)), new(new Region(0, sequence2.Length, "chr2", true)) };
        sample.EventDescs["0"] = new List<CNEventDesc>();
        sample.Kars["0"] = new Karyotype(contigs, new List<GenRange>(), _genRef.Centromeres, SexEnum.Male);
        // Apply an internal inversion
        var rnd = new Random(0);
        int contigID = 0;
        var eventData = new InternalEventData(rnd, eventPars[0], contigID, sequence1.Count());
        var eventDesc = new CNEventDesc(CNEventType.InternalInversion, 1, eventData.ToString());
        sample.EventDescs["0"].Add(eventDesc);
        sample.Kars["0"].ApplyInternalInversion(contigID, 4, 8);

        contigID = 1;
        var snvEventPars = new CNEventPars(CNEventType.SNV, 1);
        var snvEventData = new PointMutationData(rnd, snvEventPars, contigID, sequence2.Count());
        var snvEventDesc = new CNEventDesc(CNEventType.SNV, 2, snvEventData.ToString());
        sample.EventDescs["0"].Append(snvEventDesc);
        sample.Kars["0"].ApplySNV(contigID, 10, Nucleotide.A);

        var samples = new List<Sample>
            { sample };
        if (samples.Any(s => s.EventDescs.Any()))
        {
            files.WriteFasta(_genRef, samples);
        }
    }
}
