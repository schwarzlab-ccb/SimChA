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
    
    private static List<Signature> MakeSigs(List<CNEventPars> eventPs) 
        => new() { new ("TestSig", 1, eventPs) };
    
    private static List<CTreeNode> EmptyTree(CTreeNode node) => new() { node };

    [SetUp]
    public void Setup()
    {
        _rnd = new Random(0);
        _genRef = FileIO.GetGenRef("./../../../../data/hg19");
    }

    
    private Simulator GetSimulator(Type type, SimParams? simParams = null, FitParams? fitParams = null, MHParams? mhParams = null, SAParams? saParams = null)
    {
        simParams ??= new SimParams();
        fitParams ??= new FitParams();
        mhParams ??= new MHParams();
        saParams ??= new SAParams();
        return type switch
        {
            not null when type == typeof(Simulator) => new Simulator(_rnd, _genRef, simParams, fitParams),
            not null when type == typeof(MHSimulator) => new MHSimulator(_rnd, _genRef, simParams, fitParams, mhParams),
            not null when type == typeof(SASimulator) => new SASimulator(_rnd, _genRef, simParams, fitParams, saParams),
            _ => throw new ArgumentException("Unknown simulator type")
        };
    }
    
    [TestCase(typeof(Simulator)), TestCase(typeof(MHSimulator)), TestCase(typeof(SASimulator))]
    public void TestSimulatorsAll(Type simulatorType)
    {
        var sim = GetSimulator(simulatorType);
        var eventPs = new List<CNEventPars> { new(CNEventType.Pass, 1) };
        var node = new CTreeNode("root", "root", 1, 1);
        var res = sim.Simulate(node, EmptyTree(node), MakeSigs(eventPs)); 
        Assert.AreEqual(1, res.Count);
        Assert.AreEqual(1, res[0].Events.Count);
        Assert.AreEqual(46, res[0].Karyotype.CountContigs());
    }
    
    [Test]
    public void TestEmptySimulator()
    {
        
    }
}
