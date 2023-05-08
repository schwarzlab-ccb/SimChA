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
    private MCParams _mcParams;
    private Simulator _sim;
    private List<CNEventP> _eventPs;
    private FitnessParams _fitness;
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
        _fitness = new FitnessParams(1, 1, 1);
        _eventPs = new List<CNEventP> {new(CNEventType.ChromDuplication, .4), new(CNEventType.ChromDeletion, .6)};
        _mcParams = new MCParams(0, 0, 1.0, 1.0, 0.0);
        _sexXX = true;
        _geneLists = new Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>>
        {
            {GeneListType.Essentiality, new Dictionary<ChrNo, List<Gene>>()},
            {GeneListType.TumorSuppressor, new Dictionary<ChrNo, List<Gene>>()},
            {GeneListType.Oncogene, new Dictionary<ChrNo, List<Gene>>()}
        };
        _sim = new Simulator(_rnd, _fitness, _mcParams, _geneLists);
        _fitnessDict = new Dictionary<string, double> { {"0", 1.0} };
        var kar = new Karyotype(_sexXX);
        _clone = new Clone(0, -1, "0", 0, kar, 1);
    }
    
    // Taken the MakeGene method from TestFitness
    private static Gene MakeGene(ChrNo chrNo, double deltaFitness)
        => new($"G{chrNo}", new Region(0, 50, new ChrID(chrNo, false)), deltaFitness);
    
    [Test]
    public void TestPotential()
    {
        var events = new List<BaseEventData>();
        bool threshold = false;
        var listGenes = Enum.GetValues(typeof(GeneListType)).Cast<GeneListType>().ToDictionary(
            t => t,
            _ => Enum.GetValues(typeof(ChrNo)).Cast<ChrNo>().ToDictionary(chrNo => chrNo, _ => new List<Gene>()));

        listGenes[GeneListType.Oncogene][ChrNo.chr1].Add(MakeGene(ChrNo.chr1, 0.001));
        _sim = new Simulator(_rnd, _fitness, _mcParams, listGenes);
        double potential = _sim.Potential(_clone, _fitnessDict, events, ref threshold);
        Assert.AreEqual(potential,1.0,EPSILON);
    }
    
    // [Test]
    // public void TestInitEventPs()
    // {
    //     // Test 0 generation
    //     int nMutations = 0;
    //     var emptyEventPs = _sim.InitEventPs(_clone, nMutations).ToList();
    //     Assert.AreEqual(emptyEventPs.Count, nMutations);
    //     
    //     // Test that it generates the correct number of mutations needed
    //     nMutations = 5;
    //     var eventPs = _sim.InitEventPs(_clone, nMutations).ToList();
    //     Assert.AreEqual(eventPs.Count, nMutations);
    //     
    //     // Test that it will never sample 0 prob events
    //     var events = new List<CNEventP> {new(CNEventType.ChromDuplication, 1), new(CNEventType.ChromDeletion, 0)};
    //     var sigs = new List<Signature> { new("test", 1, events)};
    //     _sim = new Simulator(_rnd, _fitness, sigs, _mcParams, _geneLists);
    //     nMutations = 5;
    //     eventPs = _sim.InitEventPs(_clone, nMutations).ToList();
    //     foreach (var e in eventPs)
    //     {
    //         Assert.AreEqual(e.Type, CNEventType.ChromDuplication);
    //     }
    //     Assert.AreEqual(eventPs.Count, nMutations);
    // }

    [Test]
    public void TestInitEvents()
    {
        // Test that it will never sample 0 prob signatures
        // Init a new Simulator instance
        _sim = new Simulator(_rnd, _fitness,  _mcParams, _geneLists);
        int nMutations = 5;
        var eventData = _sim.InitEvents(_clone.Karyotype, nMutations, _eventPs);
        foreach (var data in eventData)
        {
            Assert.True(data.EventType is CNEventType.ChromDeletion or CNEventType.ChromDuplication);
        }
        Assert.AreEqual(eventData.Count, nMutations);
    }
}