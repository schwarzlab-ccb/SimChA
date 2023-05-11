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
    private List<CNEventPars> _eventPs;
    private FitnessParams _fitness;
    private Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> _geneLists;
    private Karyotype _kar;
    private const double EPSILON = 0.0000000001;
    
    [SetUp]
    public void Setup()
    {
        _rnd = new Random(0);
        _fitness = new FitnessParams(1, 1, 1);
        _eventPs = new List<CNEventPars> {new(CNEventType.ChromDuplication, .4), new(CNEventType.ChromDeletion, .6)};
        _mcParams = new MCParams(0, 0, 1.0, 1.0, 0.0);
        _geneLists = new Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>>
        {
            {GeneListType.Essentiality, new Dictionary<ChrNo, List<Gene>>()},
            {GeneListType.TumorSuppressor, new Dictionary<ChrNo, List<Gene>>()},
            {GeneListType.Oncogene, new Dictionary<ChrNo, List<Gene>>()}
        };
        _kar = new Karyotype(true);
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
        var sim = new Simulator(_rnd, _fitness, listGenes);
        double potential = sim.Potential(_mcParams, _kar, 1, events).potential;
        Assert.AreEqual(potential,0.0,EPSILON);
    }

    [Test]
    public void TestInitEvents()
    {
        // Test that it will never sample 0 prob signatures
        // Init a new Simulator instance
        var sim = new Simulator(_rnd, _fitness, _geneLists);
        const int nMutations = 5;
        var eventData = sim.InitEvents(_kar, nMutations, _eventPs);
        foreach (var data in eventData)
        {
            Assert.True(data.EventType is CNEventType.ChromDeletion or CNEventType.ChromDuplication);
        }
        Assert.AreEqual(eventData.Count, nMutations);
    }
}