using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using NUnit.Framework;
using SimChA.Data;
using SimChA.IO;
using JsonSerializer = System.Text.Json.JsonSerializer;

namespace Tests;

[TestFixture]
public class TestParsing
{
    public const string DATA_PATH = "./../../../../data/";
    public const string HG_19_PATH = DATA_PATH + "./hg19";
    public const string HG_38_PATH = DATA_PATH + "./hg38";
    private GenRef _genRef;

    [SetUp]
    public void Setup()
    {
        _genRef = FileIO.ReadGenRef(HG_19_PATH);
    }
    
    [Test]
    public void TestConfigSerialization()
    {
        var config = new SimChAConfig(new SimParams(), new ChAParams(), new FitParams());
        var options = new JsonSerializerOptions { WriteIndented = true };
        string serialized = JsonSerializer.Serialize(config, options);
        var deserialized = JsonSerializer.Deserialize<SimChAConfig>(serialized);
        Assert.NotNull(deserialized);
        // Assure that the deserialized object is the same as the original, including the nested objects
        Assert.AreEqual(config.SimParams, deserialized?.SimParams);
        Assert.AreEqual(config.Signatures, deserialized?.Signatures);
    }
    
    [Test]
    public void TestParseGeneLists()
    {
        var tsgList = _genRef.AllChrs.ToDictionary(t => t, t => new List<Gene>());

        string genesTSG = "chr1\t1\t50\tTSG1\t0.001";

        var gene1 = new Gene("TSG1", new GenRange(0, 50, "chr1"), 0.001);
        tsgList["chr1"].Add(gene1);
        var listFromString = Parsers.ParseGeneList(new StringReader(genesTSG), _genRef.AllChrs);
        Assert.AreEqual(tsgList["chr1"], listFromString);

        genesTSG += "\nchr2\t100\t5000\tTSG2\t0.01";
        var gene2 = new Gene("TSG2", new GenRange(99, 5000, "chr2"), 0.01);
        tsgList["chr2"].Add(gene2);
        listFromString = Parsers.ParseGeneList(new StringReader(genesTSG), _genRef.AllChrs);
        Assert.AreEqual(tsgList, listFromString);
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
        Assert.AreEqual(1.728, clones[3].Fitness, double.Epsilon * 10);
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
        Assert.AreEqual(2, profiles["1"].FindChrRegions("chr1").Count()); 
        Assert.AreEqual(0, profiles["1"].FindChrRegions("chr2").Count());
        Assert.AreEqual(5, profiles["1"].FindChrRegions("chr3").Count()); 
        Assert.AreEqual(1, profiles["2"].FindChrRegions("chrX").Count()); 
        Assert.AreEqual(SexType.Male, profiles["2"].Sex);
    }

    [Test]
    public void TestParseChromFile()
    {
        Assert.AreEqual("hg19", _genRef.Name);
        Assert.AreEqual(22, _genRef.AutosomesCount);
        Assert.AreEqual(46, _genRef.ChrCount(SexType.Female, true));
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
}
