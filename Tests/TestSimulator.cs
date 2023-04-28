using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SimChA.DataTypes;
using SimChA.EventData;
using SimChA.Simulation;

namespace Tests;

[TestFixture]
public class TestSimulator
{
    private Random _rnd;
    private SimParams _simParams;
    private Simulator _sim;
    private List<Signature> _signatures;
    private FitnessParams _fitnessParams;
    private Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> _geneLists;
    private Dictionary<string, double> _fitnessDict;
    
    [SetUp]
    public void Setup()
    {
        int seed = 0;
        _rnd = new Random(seed);
        _fitnessParams = new FitnessParams(1, 1, 1);
        var events = new List<CNEventP> {new(CNEventType.ChromDuplication, 1), new(CNEventType.ChromDeletion, 1)};
        _signatures = new List<Signature> { new("test", 1, events)};
        _simParams = new SimParams(
            seed, 
            true, 
            1, 
            Distribution.Uniform, 
            GenomeAssembly.hg38, 
            _fitnessParams, 
            _signatures, 
            null);
        _geneLists = new Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>>()
        {
            {GeneListType.Essentiality, new Dictionary<ChrNo, List<Gene>>()},
            {GeneListType.TumorSuppressor, new Dictionary<ChrNo, List<Gene>>()},
            {GeneListType.Oncogene, new Dictionary<ChrNo, List<Gene>>()}
        };
        _sim = new Simulator(_rnd, _simParams, _geneLists);
        _fitnessDict = new Dictionary<string, double>() { {"0", 0.0} };
    }

    [Test]
    public void TestPotential()
    {
        // TODO: implement the test for potential
    }

    [Test]
    public void TestInitEventPs()
    {
        // Test 0 generation
        int nMutations = 0;
        List<CNEventP> eventPs = _sim.InitEventPs(_signatures[0], nMutations).ToList();
        
        Assert.AreEqual(eventPs.Count(), nMutations);
        // Test that it generates the correct number of mutations needed
        var events = new List<CNEventP> {new(CNEventType.ChromDuplication, 1)};
        var sig = new Signature("test", 1, events);
        nMutations = 5;
        eventPs = _sim.InitEventPs(sig, nMutations).ToList();
        foreach (var e in eventPs)
        {
            Assert.AreEqual(e.Type, "ChromDuplication");
        }
        Assert.AreEqual(eventPs.Count(), nMutations);
        
        // Test that it will never sample 0 prob events
        events = new List<CNEventP> {new(CNEventType.ChromDuplication, 1), new(CNEventType.ChromDeletion, 0)};
        sig = new Signature("test", 1, events);
        nMutations = 5;
        eventPs = _sim.InitEventPs(sig, nMutations).ToList();
        foreach (var e in eventPs)
        {
            Assert.AreEqual(e.Type, "ChromDuplication");
        }
        Assert.AreEqual(eventPs.Count(), nMutations);
    }

    [Test]
    public void TestInitEvents()
    {
        // Test that it will never sample 0 prob signatures
        var events = new List<CNEventP> {new(CNEventType.ChromDuplication, 1)};
        var sig1 = new Signature("test1", 1, events);
        events = new List<CNEventP> {new(CNEventType.ChromDeletion, 1)};
        var sig2 = new Signature("test2", 0, events);
        var sigs = new List<Signature>() { sig1, sig2 };
        // TODO: finish this test
    }
}