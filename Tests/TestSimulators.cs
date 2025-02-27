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
    private Random _rnd;
    private GenRef _genRef;
    
    private static List<Signature> MakeSigs(List<CNEventPars> eventPs) 
        => new() { new ("TestSig", 1, eventPs) };
    
    private static List<CTreeNode> EmptyTree(CTreeNode node) => new() { node };

    private List<Signature> _baseEvs;

    [SetUp]
    public void Setup()
    {
        _rnd = new Random(0);
        _genRef = FileIO.ReadGenRef("./../../../../data/hg19");
    }
    
    private Simulator GetSimulator(Type type, SimParams? simParams = null, FitParams? fitParams = null, MHParams? mhParams = null, EvoParams? saParams = null)
    {
        simParams ??= new SimParams();
        fitParams ??= new FitParams();
        mhParams ??= new MHParams();
        saParams ??= new EvoParams();
        return type switch
        {
            not null when type == typeof(Simulator) => new Simulator(_rnd, _genRef, simParams, fitParams),
            not null when type == typeof(MHSimulator) => new MHSimulator(_rnd, _genRef, simParams, fitParams, mhParams),
            not null when type == typeof(EvoSimulator) => new EvoSimulator(_rnd, _genRef, simParams, fitParams, saParams),
            _ => throw new ArgumentException("Unknown simulator type")
        };
    }
    
    [TestCase(typeof(Simulator)), TestCase(typeof(MHSimulator)), TestCase(typeof(EvoSimulator))]
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
    
    [TestCase(typeof(Simulator)), TestCase(typeof(MHSimulator)), TestCase(typeof(EvoSimulator))]
    public void TestEmptySimulator(Type simulatorType)
    {
        var sim = GetSimulator(simulatorType);
        var eventPs = new List<CNEventPars> { new(CNEventType.ChromDeletion, 1) };
        int dist = 50;
        var node = new CTreeNode("root", "root", dist, 1);
        var res = sim.Simulate(node, EmptyTree(node), MakeSigs(eventPs)); 
        Assert.AreEqual(1, res.Count);
        Assert.AreEqual(dist, res[0].Events.Count);
    }
    
    [TestCase(typeof(Simulator)), TestCase(typeof(MHSimulator)), TestCase(typeof(EvoSimulator))]
    public void TestSimulatorsBase(Type simulatorType)
    {
        var sim = GetSimulator(simulatorType);
        List<CNEventPars> eventPs = new()
        {
            new CNEventPars(CNEventType.ChromDuplication, 1),
            new CNEventPars(CNEventType.ChromDeletion, 1),
            new CNEventPars(CNEventType.ArmDeletion, 1),
            new CNEventPars(CNEventType.ArmDuplication, 1),
            new CNEventPars(CNEventType.InternalDeletion, 1, .1),
            new CNEventPars(CNEventType.InternalDuplication, 1, .1),
            new CNEventPars(CNEventType.InternalInversion, 1, .1),
            new CNEventPars(CNEventType.InvertedDuplication, 1, .1),
            new CNEventPars(CNEventType.TailDeletion, 1, .3),
            new CNEventPars(CNEventType.TailDuplication, 1, .3),
            new CNEventPars(CNEventType.CentromereBoundDeletion, 1, .2),
            new CNEventPars(CNEventType.CentromereBoundDuplication, 1, .2),
        };
        int dist = 100;
        var root = new CTreeNode("root", "root", 0, 1);
        var tree = Enumerable.Range(1, 10)
            .Select(i => new CTreeNode("clone" + i, "root", dist, 1));
        var res = sim.Simulate(root, tree.ToList(), MakeSigs(eventPs)); 
        Assert.AreEqual(10, res.Count);
        Assert.AreEqual(dist, res[0].Events.Count);
    }
}