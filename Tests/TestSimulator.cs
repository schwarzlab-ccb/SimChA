using System;
using System.Collections.Generic;
using NUnit.Framework;
using SimChA.Data;
using SimChA.EventData;
using SimChA.IO;
using SimChA.Simulation;

namespace Tests;

[TestFixture]
public class TestSimulator
{
    private Random _rnd;
    private SimParams _simParams;
    private MHParams _mhParams;
    private List<CNEventPars> _eventPs;
    private FitParams _fitParams;
    private GenRef _genRef;
    private Karyotype _kar;
    private const double EPSILON = 0.0000000001;    
    
    [SetUp]
    public void Setup()
    {
        _rnd = new Random(0);
        _simParams = new SimParams();
        _fitParams = new FitParams(0.9, 0.05, 2);
        _eventPs = new List<CNEventPars> {new(CNEventType.ChromDuplication, .4), new(CNEventType.ChromDeletion, .6)};
        _genRef = FileIO.GetGenRef("./../../../../data/hg19");
        _mhParams = new MHParams(0, 0, 1.0, true, 1.0, 0.0);
        _kar = new Karyotype(_genRef, SexType.Female);
    }
}
