// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

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

namespace Tests;

[TestFixture]
public class TestIO
{
    [Test]
    public void TestConfigSerialization()
    {
        var fit = new FitnessParams(0.001f, 0.01f, 0.000_1f);
        var simParams = new SimParams(0, true, 1, Distribution.Uniform, GenomeAssembly.hg38, fit, null, null);
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
        Assert.AreEqual(CNEventType.WholeGenomeDoubling, res.Signatures!.First().Events!.First().Type);
        res = Parsers.ParseSimParams(@"{""Signatures"": [{""Name"": ""test"", ""Prob"": 1, ""Events"": [{""Type"": ""InternalInversion"", ""Prob"": 0.1, ""Params"": {""Mean"": 0.1}}]}]}");
        Assert.AreEqual(0.1, res.Signatures!.First().Events!.First().Params!["Mean"], 0.000001);
        res = Parsers.ParseSimParams(@"{""Assembly"":""hg38""}");
        Assert.AreEqual(GenomeAssembly.hg38, res.Assembly);
    }

    [Test]
    public void TestNewickParser()
    {
        string[] newickTestStrings= {
            "",
            "(,,(,));",
            "(A,B,(C,D));",
            "(A,B,(C,D)E)F;",
            "(:1,:2,(:3,:4):5);",
            "(:1,:2,(:3,:4):5):0;",
            "(A:1,B:2,(C:3,D:4):5);",
            "(A:1,B:2,(C:3,D:4)E:5)F;",
            "((B:2,(C:3,D:4)E:5)F:1)A;",
            "(((RRib7Met_A24D-0021_CRUK_PC_0021_M3:8.00000,(RDiaphragmaticMet_A24C-0021_CRUK_PC_0021_M2:4.00000,XiphoidMet_A24E-0021_CRUK_PC_0021_M4:1.00000)internal_2:1.00000)"+
            "internal_3:4.00000,RAxillaryLNMet_A24A-0021_CRUK_PC_0021_M1:30.00000)internal_1:44.00000,diploid:0.00000):0.00000;",
            "(((6506-1:12, ((5445-1:2)3158-1:6, 6419-1:21)1094-0:6)400-0:13, (6330-1:26, ((((3448-1:6, (5382-1:13, 6238-1:4)1758-0:1)1725-0:13, "+
            "((4011-1:15, (((5405-1:2)4383-1:4, (6001-1:8, ((6681-1:2)4302-1:3, 6047-1:5)3776-0:6)2827-0:2)2475-0:3, (6132-1:1, (6475-1:1, 6626-1:3)"+
            "5398-0:1)5135-0:6)1706-0:8)744-0:4, (5183-1:2, 5279-1:4)3794-0:7)564-0:5)419-0:9, ((6920-1:15, ((((5568-1:12, ((6921-1:6)4427-1:11, (6310-1:1)"+
            "5786-1:17)1359-0:3)1142-0:1, (((6838-1:1)6648-1:6, ((5473-1:5, 6034-1:3)4325-0:2, (5994-1:4, 6142-1:3)4129-0:3)3514-0:1)3213-0:8, ((6592-1:4, "+
            "(5555-1:2, (6453-1:3, 6654-1:3)4882-0:1)4282-0:1, (((6846-1:1)6421-1:1)5865-1:3, (6799-1:1, 6905-1:1)5882-0:7)4520-0:2)4045-0:5, (6465-1:9, "+
            "((6509-1:1)5645-1:4, (6869-1:1)6735-1:2, (5995-1:1, (6790-1:2, 6854-1:2)6342-1:2, 6910-1:2)5494-0:4)3454-0:3)3006-0:1)2107-0:10)1238-0:2)1134-0:1, "+
            "((4565-1:1, 6427-1:3, (6168-1:2, 6445-1:2)5337-0:3)3756-0:9, (5161-1:4, (5519-1:5, 5879-1:6)3831-0:2)2879-0:10)1329-0:2)1057-0:2, (6781-1:13, "+
            "(6650-1:8, ((6230-1:1)5616-1:1, 6617-1:5)5139-0:3, (6169-1:6, (5806-1:1, 6146-1:1, 6392-1:1)5754-0:5)3918-0:1, ((6173-1:1)6154-1:1, 6746-1:2)"+
            "5360-0:2)3539-0:11)1082-0:3)698-0:1)636-0:3, ((6302-1:8, (5831-1:4, (6301-1:1, 6643-1:3)4570-0:1)4183-0:4)2190-0:7, ((5897-1:2, 6883-1:7)4425-0:7, "+
            "(6614-1:2, 6702-1:3)4081-0:10)1627-0:6)623-0:2, (6398-1:13, 6824-1:17)878-0:6)598-0:8)178-0:3, ((5282-1:12, (6018-1:8, 6904-1:8)1923-0:4)1160-0:2, "+
            "(6582-1:14, (5385-1:1, 6129-1:1, 6745-1:5)4070-0:5)1882-0:2)858-0:9)130-0:3)55-0:5, ((6794-1:1)5640-1:22, ((5707-1:17, ((4555-1:4, 6122-1:6, 6676-1:8)3281-0:2, "+
            "((5051-1:4, 6459-1:11)3563-0:2, ((6841-1:1)5163-1:1, 6756-1:2)4571-0:9)2687-0:4)2183-0:7)635-0:2, (6740-1:6, 6830-1:3)3456-0:12)571-0:1)535-0:7)1-0:1, 0-0:0):0;"
        };
        
        int[] cloneCounts = {
            0,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            9,
            174
        };

        for (int i = 0; i < newickTestStrings.Length; i++)
        {
            var clones = Parsers.ParseNewick(newickTestStrings[i], true);
            foreach(var clone in clones)
            {
                Console.WriteLine($"String {i}: " +
                                  $"ID={clone.CloneId}, " +
                                  $"Name={clone.Name}, " +
                                  $"ParentID={clone.ParentId}, " +
                                  $"Mutations={clone.DistToParent}, " +
                                  $"ChildrenCount={clone.ChildrenIDs.Count}");
                foreach(int childrenID in clone.ChildrenIDs)
                {
                    Console.WriteLine($"\tChildID:{childrenID}");
                }
                Console.WriteLine();
            }
            Assert.AreEqual(clones.Count, cloneCounts[i]);
            Console.WriteLine();
        }
    }

    [Test]
    public void TestReadGeneLists()
    {
        var tsgList = Enum.GetValues<ChrNo>().ToDictionary(t => t, t => new List<Gene>());
        
        string genesTSG = "chr1\t1\t50\tTSG1\t0.001";
        
        var gene1 = new Gene("TSG1", new GenRange(0, 50, ChrNo.chr1), 0.001);
        tsgList[ChrNo.chr1].Add(gene1);
        var listFromString = Parsers.ParseGeneList(new StringReader(genesTSG), true);
        Assert.AreEqual(tsgList, listFromString);
        
        genesTSG += "\nchr2\t100\t5000\tTSG2\t0.01";
        var gene2 = new Gene("TSG2", new GenRange(99, 5000, ChrNo.chr2), 0.01);
        tsgList[ChrNo.chr2].Add(gene2);
        listFromString = Parsers.ParseGeneList(new StringReader(genesTSG), true);
        Assert.AreEqual(tsgList, listFromString);
    }

    [Test]
    public void TestWrite()
    {
        var projectPath = Path.GetDirectoryName(Path.GetDirectoryName(Path.GetDirectoryName(TestContext.CurrentContext.TestDirectory)));
        var files = new FileIO(projectPath + "/out");
        var kar = new Karyotype(false);
        var ceParams = new Dictionary<string, double> {{"Size", 2000000}};
        string eventDesc = kar.ApplyCNEvent(new Random(48), new CNEventP(CNEventType.Rigma, 1.0, ceParams));
        Console.WriteLine(eventDesc);
        var clone = new Clone(1, -1, "test", 0, kar, 1);
        files.WriteClones(new List<Clone> {clone});
    }
}