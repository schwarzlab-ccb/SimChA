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
    private MCParams _mcParams;
    private Simulator _sim;
    private List<Signature> _signatures;
    private FitnessParams _fitnessParams;
    private Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> _geneLists;
    private Dictionary<string, double> _fitnessDict;
    private int _seed;
    private bool _sexXX;
    private Clone _clone;
    private const double EPSILON = 0.0000000001;
    
    [SetUp]
    public void Setup()
    {
        _seed = 0;
        _rnd = new Random(_seed);
        _fitnessParams = new FitnessParams(1, 1, 1);
        var events = new List<CNEventP> {new(CNEventType.ChromDuplication, 1), new(CNEventType.ChromDeletion, 1)};
        _signatures = new List<Signature> { new("test", 1, events)};
        _mcParams = new MCParams(0, 0, 1.0, 1.0, 0.0);
        _sexXX = true;
        _simParams = new SimParams(
            _seed, 
            _sexXX, 
            1, 
            Distribution.Uniform, 
            GenomeAssembly.hg38, 
            _fitnessParams, 
            _signatures, 
            _mcParams);
        _geneLists = new Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>>()
        {
            {GeneListType.Essentiality, new Dictionary<ChrNo, List<Gene>>()},
            {GeneListType.TumorSuppressor, new Dictionary<ChrNo, List<Gene>>()},
            {GeneListType.Oncogene, new Dictionary<ChrNo, List<Gene>>()}
        };
        _sim = new Simulator(_rnd, _simParams, _geneLists);
        _fitnessDict = new Dictionary<string, double>() { {"0", 1.0} };
        var kar = new Karyotype(_sexXX);
        _clone = new Clone(1, -1, "test", 0, kar, 1);
    }
    // Taken the MakeGene method from TestFitness
    private static Gene MakeGene(ChrNo chrNo, double deltaFitness)
        => new($"G{chrNo}", new Region(0, 50, new ChrID(chrNo, false)), deltaFitness);
    [Test]
    public void TestPotential()
    {
        List<BaseEventData> events = new List<BaseEventData>();
        bool threshold = false;
        var listGenes = Enum.GetValues(typeof(GeneListType)).Cast<GeneListType>().ToDictionary(
            t => t,
            _ => Enum.GetValues(typeof(ChrNo)).Cast<ChrNo>().ToDictionary(chrNo => chrNo, _ => new List<Gene>()));

        listGenes[GeneListType.Oncogene][ChrNo.chr1].Add(MakeGene(ChrNo.chr1, 0.001));
        _sim = new Simulator(_rnd, _simParams, listGenes);
        double potential = _sim.Potential(_clone, _fitnessDict, events, ref threshold);
        Assert.AreEqual(potential,1.0,EPSILON);
    }
    
    [Test]
    public void TestInitEventPs()
    {
        // Test 0 generation
        int nMutations = 0;
        var emptyEventPs = _sim.InitEventPs(nMutations).ToList();
        Assert.AreEqual(emptyEventPs.Count(), nMutations);
        
        // Test that it generates the correct number of mutations needed
        nMutations = 5;
        var eventPs = _sim.InitEventPs(nMutations).ToList();
        Assert.AreEqual(eventPs.Count(), nMutations);
        
        // Test that it will never sample 0 prob events
        var events = new List<CNEventP> {new(CNEventType.ChromDuplication, 1), new(CNEventType.ChromDeletion, 0)};
        var sigs = new List<Signature> { new("test", 1, events)};
        var simParams = new SimParams(
            _seed, 
            true, 
            1, 
            Distribution.Uniform, 
            GenomeAssembly.hg38, 
            _fitnessParams, 
            sigs, 
            _mcParams);
        _sim = new Simulator(_rnd, simParams, _geneLists);
        nMutations = 5;
        eventPs = _sim.InitEventPs(nMutations).ToList();
        foreach (var e in eventPs)
        {
            Assert.AreEqual(e.Type, CNEventType.ChromDuplication);
        }
        Assert.AreEqual(eventPs.Count(), nMutations);
        Assert.AreEqual(_sim.SelectedSignatures.Count(), nMutations);
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
        // Init a new Simulator instance
        var sex = true;
        var simParams = new SimParams(
            _seed, 
            sex, 
            1, 
            Distribution.Uniform, 
            GenomeAssembly.hg38, 
            _fitnessParams, 
            sigs, 
            _mcParams);
        _sim = new Simulator(_rnd, simParams, _geneLists);
        int nMutations = 5;
        var eventData = _sim.InitEvents(_clone, nMutations);
        foreach (var data in eventData)
        {
            Assert.AreEqual(data.EventType, CNEventType.ChromDuplication);
        }
        Assert.AreEqual(eventData.Count(), nMutations);
    }
}