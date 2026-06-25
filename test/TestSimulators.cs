using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SimChA.Data;
using SimChA.EventData;
using SimChA.IO;
using SimChA.Simulation;

namespace Tests;

[TestFixture]
public class TestSimulators
{
    private Random _rnd = null!;
    private RefGen _refGen = null!;
    
    private static List<Signature> MakeSigs(List<CNEventPars> eventPs) 
        => [new("TestSig", 1, eventPs)];
    
    private static List<CTreeNode> EmptyTree(CTreeNode node) => [node];
    
    [SetUp]
    public void Setup()
    {
        _rnd = new Random(0);
        _refGen = FileIO.ReadGenRef(TestParsing.DATA_PATH, TestParsing.HG_19, TestParsing.GENE_SET);
    }
    
    private Simulator GetSimulator(Type type, SimParams? simParams = null, FitParams? fitParams = null, EvoParams? saParams = null)
    {
        simParams ??= new SimParams();
        fitParams ??= new FitParams(1, 1, 1);
        saParams ??= new EvoParams(1, 10);
        return type switch
        {
            not null when type == typeof(Simulator) => new Simulator(_rnd, _refGen, simParams, fitParams),
            not null when type == typeof(MatchSimulator) => new MatchSimulator(_rnd, _refGen, simParams, fitParams, saParams),
            not null when type == typeof(EvoSimulator) => new EvoSimulator(_rnd, _refGen, simParams, fitParams, saParams),
            _ => throw new ArgumentException("Unknown simulator type")
        };
    }
    
    [TestCase(typeof(Simulator)), TestCase(typeof(MatchSimulator)), TestCase(typeof(EvoSimulator))]
    public void TestSimulatorsPass(Type simulatorType)
    {
        var sim = GetSimulator(simulatorType);
        var eventPs = new List<CNEventPars> { new(CNEventType.Pass, 1) };
        var node = new CTreeNode("root", "root", 1, 1);
        var res = sim.Simulate(node, EmptyTree(node), MakeSigs(eventPs)); 
        Assert.AreEqual(1, res.Count);
        Assert.AreEqual(1, res[0].Events.Count);
        Assert.AreEqual(46, res[0].Karyotype.CountContigs());
    }
    
    [TestCase(typeof(Simulator)), TestCase(typeof(MatchSimulator)), TestCase(typeof(EvoSimulator))]
    public void TestSeededSimulator(Type simulatorType)
    {
        var sim = GetSimulator(simulatorType);
        // Seed karyotype distinct from a diploid genome (46 contigs -> 92 after WGD).
        var seedKar = new Karyotype(_refGen, SexType.Male);
        seedKar.ApplyWGD();
        var eventPs = new List<CNEventPars> { new(CNEventType.Pass, 1) };
        var node = new CTreeNode("root", "root", 1, 1);
        var res = sim.Simulate(node, EmptyTree(node), MakeSigs(eventPs), seedKar);
        Assert.AreEqual(1, res.Count);
        // Pass events keep the karyotype unchanged, so it must reflect the seed, not a reset diploid.
        Assert.AreEqual(92, res[0].Karyotype.CountContigs());
    }

    [TestCase(typeof(Simulator)), TestCase(typeof(MatchSimulator)), TestCase(typeof(EvoSimulator))]
    public void TestSeededSimulatorIsCopied(Type simulatorType)
    {
        var sim = GetSimulator(simulatorType);
        var seedKar = new Karyotype(_refGen, SexType.Male);
        var eventPs = new List<CNEventPars> { new(CNEventType.ChromDeletion, 1) };
        var node = new CTreeNode("root", "root", 5, 1);
        var res = sim.Simulate(node, EmptyTree(node), MakeSigs(eventPs), seedKar);
        // The original seed must not be mutated by the simulation.
        Assert.AreEqual(46, seedKar.CountContigs());
        Assert.AreEqual(1, res.Count);
    }

    [TestCase(typeof(Simulator)), TestCase(typeof(MatchSimulator)), TestCase(typeof(EvoSimulator))]
    public void TestEmptySimulator(Type simulatorType)
    {
        var sim = GetSimulator(simulatorType);
        var eventPs = new List<CNEventPars> { new(CNEventType.ChromDeletion, 1) };
        const int dist = 50;
        var node = new CTreeNode("root", "root", dist, 1);
        var res = sim.Simulate(node, EmptyTree(node), MakeSigs(eventPs)); 
        Assert.AreEqual(1, res.Count);
        Assert.AreEqual(dist, res[0].Events.Count);
    }
    
    [TestCase(typeof(Simulator)), TestCase(typeof(MatchSimulator)), TestCase(typeof(EvoSimulator))]
    public void TestSimulatorsBase(Type simulatorType)
    {
        var sim = GetSimulator(simulatorType);
        List<CNEventPars> eventPs =
        [
            new(CNEventType.ChromDuplication, 1),
            new(CNEventType.ChromDeletion, 1),
            new(CNEventType.ArmDeletion, 1),
            new(CNEventType.ArmDuplication, 1),
            new(CNEventType.InternalDeletion, 1, .1),
            new(CNEventType.InternalDuplication, 1, .1),
            new(CNEventType.InternalInversion, 1, .1),
            new(CNEventType.InvertedDuplication, 1, .1),
            new(CNEventType.TailDeletion, 1, .3),
            new(CNEventType.TailDuplication, 1, .3),
            new(CNEventType.CentromereBoundDeletion, 1, .2),
            new(CNEventType.CentromereBoundDuplication, 1, .2)
        ];
        int dist = 10;
        var root = new CTreeNode("root", "root", 0, 1);
        var tree = Enumerable.Range(1, 2)
            .Select(i => new CTreeNode("clone" + i, "root", dist, 1)).ToList();
        var res = sim.Simulate(root, tree, MakeSigs(eventPs)); 
        Assert.AreEqual(tree.Count, res.Count);
        Assert.AreEqual(dist, res[0].Events.Count);
    }
}