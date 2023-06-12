using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SimChA.DataTypes;
using SimChA.EventData;
using SimChA.Simulation;
using SimChA.Computation;

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
        _eventPs = new List<CNEventPars> {new(CNEventType.ChromDuplication, .5), new(CNEventType.ChromDeletion, .5)};
        _mcParams = new MCParams(0, 0, 1.0, 1.0, 0.0);
        _geneLists = Enum.GetValues(typeof(GeneListType)).Cast<GeneListType>().ToDictionary(
            gl => gl, gl => Enum.GetValues(typeof(ChrNo)).Cast<ChrNo>().ToDictionary(
            c => c, c => new List<Gene> {MakeGene(c)}));
        Fitness.SetStartingParams(_geneLists, new FitnessParams(0.1, 0.1, 0.1));
        HGRef.Assembly = GenomeAssembly.hg19;
        _kar = new Karyotype(false);
    }
    
    // Taken the MakeGene method from TestFitness
    private Gene MakeGene(ChrNo chrNo)
        => new($"G{chrNo}", new Region(0, 50, new ChrID(chrNo, false)),  0.0);
    
    [Test]
    public void TestPotential()
    {
        var events = new List<BaseEventData>();
        var sim = new MCSimulator(_rnd, _geneLists, _mcParams);
        double potential = sim.Potential(_kar, 1, events).potential;
        Assert.AreEqual(0.0, potential,EPSILON);
    }

    [Test]
    public void TestInitEvents()
    {
        var sim = new Simulator(_rnd, _geneLists);
        const int nMutations = 5;
        var eventData = sim.InitEvents(_kar, nMutations, _eventPs);
        foreach (var data in eventData)
        {
            Assert.True(data.EventType is CNEventType.ChromDeletion or CNEventType.ChromDuplication);
        }
        Assert.AreEqual(eventData.Count, nMutations);
    }
}