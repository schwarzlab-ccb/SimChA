using System;
using System.Collections.Generic;
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
    
    private CTreeNode _baseNode;

    private static List<Signature> MakeSigs(List<CNEventPars> eventPs) 
        => new() { new ("TestSig", 1, eventPs) };
    
    private static List<CTreeNode> EmptyTree(CTreeNode node) => new() { node };

    [SetUp]
    public void Setup()
    {
        _rnd = new Random(0);
        _genRef = FileIO.GetGenRef("./../../../../data/hg19");
        _baseNode = new CTreeNode("root", "root", 1, 1);
    }

    [Test]
    public void TestSimulator()
    {
        var sim = new Simulator(_rnd, _genRef, new SimParams(), new FitParams());
        var eventPs = new List<CNEventPars> { new(CNEventType.Pass, 1) };
        var res = sim.Simulate(_baseNode, EmptyTree(_baseNode), MakeSigs(eventPs)); 
        Assert.AreEqual(1, res.Count);
        Assert.AreEqual(1, res[0].Events.Count);
        Assert.AreEqual(46, res[0].Karyotype.CountContigs());
    }

    [Test]
    public void TestMHSimulator()
    {
        var sim = new MHSimulator(_rnd, _genRef, new SimParams(), new FitParams(),  new MHParams());
        var eventPs = new List<CNEventPars> { new(CNEventType.Pass, 1) };
        var res = sim.Simulate(_baseNode, EmptyTree(_baseNode), MakeSigs(eventPs)); 
        Assert.AreEqual(res.Count, 1);
        Assert.AreEqual(res[0].Events.Count, 1);
        Assert.AreEqual(res[0].Karyotype.CountContigs(), 46);
    }
    
    [Test]
    public void TestSASimulator()
    {
        var sim = new SASimulator(_rnd, _genRef, new SimParams(), new FitParams(), new SAParams());
        var eventPs = new List<CNEventPars> { new(CNEventType.Pass, 1) };
        var res = sim.Simulate(_baseNode, EmptyTree(_baseNode), MakeSigs(eventPs)); 
        Assert.AreEqual(res.Count, 1);
        Assert.AreEqual(res[0].Events.Count, 1);
        Assert.AreEqual(res[0].Karyotype.CountContigs(), 46);
    }
}
