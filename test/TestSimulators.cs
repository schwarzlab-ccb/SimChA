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
    public void TestMaxWgdLimit(Type simulatorType)
    {
        const int maxWgd = 2;
        var simParams = new SimParams(MaxWGD: maxWgd);
        var sim = GetSimulator(simulatorType, simParams);
        // A mix where WGD is possible but not guaranteed, so most samples already satisfy the
        // limit and the few that exceed it are restarted until they do.
        List<CNEventPars> eventPs =
        [
            new(CNEventType.WholeGenomeDoubling, 1),
            new(CNEventType.Pass, 4)
        ];
        var node = new CTreeNode("root", "root", 10, 1);
        var res = sim.Simulate(node, EmptyTree(node), MakeSigs(eventPs));
        int wgdCount = res[0].Events.Count(e => e.EventData.EventType == CNEventType.WholeGenomeDoubling);
        Assert.AreEqual(1, res.Count);
        Assert.AreEqual(10, res[0].Events.Count);
        Assert.LessOrEqual(wgdCount, maxWgd);
    }

    [Test]
    public void TestMaxWgdDisabledByDefault()
    {
        // Default SimParams.MaxWGD is -1 (no limit), so a WGD-only signature is never restarted.
        // The basic Simulator applies every drawn event, producing exactly `dist` WGDs.
        var sim = GetSimulator(typeof(Simulator));
        var eventPs = new List<CNEventPars> { new(CNEventType.WholeGenomeDoubling, 1) };
        var node = new CTreeNode("root", "root", 3, 1);
        var res = sim.Simulate(node, EmptyTree(node), MakeSigs(eventPs));
        int wgdCount = res[0].Events.Count(e => e.EventData.EventType == CNEventType.WholeGenomeDoubling);
        Assert.AreEqual(3, wgdCount);
    }

    [Test]
    public void TestMaxWgdUnsatisfiableThrows()
    {
        // The basic Simulator applies every drawn event, so a WGD-only signature always produces
        // `dist` WGDs. A limit of 0 is then impossible to satisfy and the run aborts rather than
        // looping forever.
        var simParams = new SimParams(MaxWGD: 0);
        var sim = GetSimulator(typeof(Simulator), simParams);
        var eventPs = new List<CNEventPars> { new(CNEventType.WholeGenomeDoubling, 1) };
        var node = new CTreeNode("root", "root", 1, 1);
        Assert.Throws<Exception>(() => sim.Simulate(node, EmptyTree(node), MakeSigs(eventPs)));
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