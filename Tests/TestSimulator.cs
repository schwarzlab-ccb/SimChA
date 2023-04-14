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
    }
}