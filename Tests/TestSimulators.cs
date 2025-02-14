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
        var sim = new Simulator(_rnd, _genRef, new SimParams(), new FitParams(0.9, 0.05, 2));
        var eventPs = new List<CNEventPars> {new(CNEventType.ChromDuplication, .4), new(CNEventType.ChromDeletion, .6)};
        var res = sim.Simulate(_baseNode, EmptyTree(_baseNode), MakeSigs(eventPs)); 
    }

    [Test]
    public void TestMHSimulator()
    {
        var mhParams = new MHParams(0, 0, 1.0, true, 1.0, 0.0);
        var sim = new MHSimulator(_rnd, _genRef, new SimParams(), new FitParams(0.9, 0.05, 2), mhParams);
        var eventPs = new List<CNEventPars> {new(CNEventType.ChromDuplication, .4), new(CNEventType.ChromDeletion, .6)};
        var res = sim.Simulate(_baseNode, EmptyTree(_baseNode), MakeSigs(eventPs)); 
    }
    
    [Test]
    public void TestSASimulator()
    {
        var saParams = new SAParams(1, 1, 1, 10, false);
        var sim = new SASimulator(_rnd, _genRef, new SimParams(), new FitParams(0.9, 0.05, 2), saParams);
        var eventPs = new List<CNEventPars> {new(CNEventType.ChromDuplication, .4), new(CNEventType.ChromDeletion, .6)};
        var res = sim.Simulate(_baseNode, EmptyTree(_baseNode), MakeSigs(eventPs)); 
    }
}
