// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SimChA.EventData;

namespace Tests;

public class TestEventData
{
    private Random _rnd;
    private Dictionary<string, double> _pars;
    
    [SetUp]
    public void Setup()
    {
        _rnd = new Random(42);
    }
    
    [Test]
    public void TestContigEventData()
    {
        var eventP = new CNEventPars(CNEventType.ChromDeletion, 1, 1_000_000);
        var eventData = new ContigEventData(eventP, 1);
        Assert.AreEqual("contig:1", eventData.ToString());
    }
    
    [Test]
    public void TestBaseEventData()
    {
        var eventP = new CNEventPars(CNEventType.WholeGenomeDoubling, 1, 1_000_000);
        var eventData = new BaseEventData(eventP);
        Assert.AreEqual("", eventData.ToString());
    }

    [Test]
    public void TestInternalEventData()
    {
        const long len = 1_000_000;
        var eventP = new CNEventPars(CNEventType.InternalDuplication, 1, 1_000_000);
        var eventData = new InternalEventData(_rnd, eventP, 0, len);
        Assert.GreaterOrEqual(eventData.Start, 0);
        Assert.LessOrEqual(eventData.Start, len);
        Assert.Greater(eventData.End,eventData.Start);
        Assert.LessOrEqual(eventData.Start, len);
    }

    [Test]
    public void TestTailEventData()
    {
        const long len = 1_000_000;
        var eventP = new CNEventPars(CNEventType.TailDeletion, 1, 1_000_000);
        var eventData = new InternalEventData(_rnd, eventP, 0, len);
        Assert.GreaterOrEqual(eventData.Start, 0);
        Assert.LessOrEqual(eventData.Start, len);
        Assert.Greater(eventData.End,eventData.Start);
        Assert.LessOrEqual(eventData.Start, len);
    }

    [Test]
    public void TestPairEventData()
    {
        const long lenA = 1_000_000;
        const long lenB = 10_000_000;
        var eventP = new CNEventPars(CNEventType.Translocation, 1, 1_000_000);
        var eventData = new PairEventData(_rnd, eventP, 0, lenA, 1, lenB);
        Assert.GreaterOrEqual(eventData.PosA, 0);
        Assert.LessOrEqual(eventData.PosA, lenA);
        Assert.GreaterOrEqual(eventData.PosB, 0);
        Assert.LessOrEqual(eventData.PosB, lenB);
        Assert.IsFalse(eventData.Inverted);
    }

    [Test]
    public void TestPyrgo()
    {
        const long len = 1_000_000;
        var eventP = new CNEventPars(CNEventType.Pyrgo, 1, 1_000_000, 10);
        var eventData = new PyrgoEventData(_rnd, eventP, 0, len);
        foreach (var frag in eventData.FragmentsList)
        {
            Assert.GreaterOrEqual(frag.start, 0);
            Assert.LessOrEqual(frag.start + frag.length, len);
        }
    }
    
    [Test]
    public void TestRigma()
    {
        const long len = 1_000_000;
        var eventP = new CNEventPars(CNEventType.Rigma, 1, 1_000_000, 10);
        var eventData = new RigmaEventData(_rnd, eventP, 0, len);
        Assert.GreaterOrEqual(eventData.Start, 0);
        foreach (long stop in eventData.StopsList)
        {
            Assert.LessOrEqual(stop + eventData.Start, len);
        }
    }

    [Test]
    public void TestTemplatedEvent()
    {
        var eventP = new CNEventPars(CNEventType.TIBridge, 1, 1_000_000, 10);
        var frags = new List<(int, long)>
        {
            (0, 100_000),
            (1, 1_000_000),
            (2, 2_000_000),
            (3, 3_000_000)
        };
        var eventData = new TemplatedEventData(_rnd, eventP, frags);
        foreach (var frag in eventData.Frags)
        {
            Assert.GreaterOrEqual(frag.start, 0);
            Assert.LessOrEqual(frag.start + frag.len, frags[frag.id].Item2);
        }
    }

    [Test]
    public void TestChromoplexy()
    {
        var eventP = new CNEventPars(CNEventType.Chromoplexy, 1, 1_000_000, 10);
        var frags = new List<(int, long)>
        {
            (0, 5_000_000),
            (1, 10_000_000),
            (2, 20_000_000),
            (3, 30_000_000)
        };
        var eventData = new ChromoplexyEventData(_rnd, eventP, frags);
        Assert.AreEqual(eventData.Stops.Sum(s => s.Count + 1), eventData.Sequence.Count);
        Assert.AreEqual(eventData.ContigIds.Count - 1, eventData.Breakpoints.Count);
        Console.WriteLine(eventData);
    }

    [Test]
    public void TestChromothripsis([Values] IntEdgeCases seed)
    {
        const long len = 100_000_000;
        var eventP = new CNEventPars(CNEventType.Chromothripsis, 1, 1_000_000);
        var eventData = new ChromothripsisEventData(new Random((int) seed), eventP, 0, len);
        foreach (long stop in eventData.StopsList)
        {
            Assert.Greater(stop, 0);
            Assert.Less(stop, len);
        }
        Assert.GreaterOrEqual(eventData.StopsList.Count + 1,eventData.SelectionList.Count);
    }
}