// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SimChA.Data;
using SimChA.EventData;

namespace Tests;

public class TestEventData
{
    private Random _rnd = null!;
    
    [SetUp]
    public void Setup()
    {
        _rnd = new Random(42);
    }
    
    [Test]
    public void TestContigEventData()
    {
        var eventP = new CNEventPars(CNEventType.ContigDeletion, 1);
        var eventData = new ContigEventData(eventP, 1, 100_000_000);
        Assert.AreEqual("contig:1;length:100000000;", eventData.EventDesc());
    }
    
    [Test]
    public void TestBaseEventData()
    {
        var eventP = new CNEventPars(CNEventType.WholeGenomeDoubling, 1, 0.1);
        var eventData = new BaseEventData(eventP);
        Assert.AreEqual("", eventData.EventDesc());
    }

    [Test]
    public void TestInternalEventData()
    {
        const long len = 1_000_000;
        var eventP = new CNEventPars(CNEventType.InternalDuplication, 1, 0.1);
        var eventData = new InternalEventData(_rnd, eventP, 0, len, new TerminalArm(true, len, len));
        Assert.GreaterOrEqual(eventData.Start, 0);
        Assert.LessOrEqual(eventData.Start, len);
        Assert.Greater(eventData.End,eventData.Start);
        Assert.LessOrEqual(eventData.Start, len);
    }

    [TestCase(true)]
    [TestCase(false)]
    public void TestInternalEventDrawsLengthThenPlacesUniformlyWithinSelectedArm(bool direction)
    {
        const long contigLen = 3_000_000;
        const long usableLen = 1_000_000;
        // The fitted scale reaches past the arm proper, into the middle of the centromere.
        var arm = new TerminalArm(direction, usableLen + 50_000, usableLen);
        var eventP = new CNEventPars(CNEventType.InternalDuplication, 1, 0.1);
        var normalizedStarts = new List<double>();

        for (int i = 0; i < 20000; i++)
        {
            var eventData = new InternalEventData(_rnd, eventP, 0, contigLen, arm);
            long armStart = direction ? 0 : contigLen - usableLen;
            long lastStart = armStart + usableLen - (eventData.End - eventData.Start);

            Assert.GreaterOrEqual(eventData.Start, armStart);
            Assert.LessOrEqual(eventData.End, armStart + usableLen);
            if (lastStart > armStart)
            {
                normalizedStarts.Add(
                    (eventData.Start - armStart) / (double) (lastStart - armStart));
            }
        }

        Assert.AreEqual(0.5, normalizedStarts.Average(), 0.02);
    }

    [TestCase(true)]
    [TestCase(false)]
    public void TestTelomereEventUsesFixedAlphaBetaLength(bool direction)
    {
        const long contigLen = 3_000_000;
        const long armLen = 1_000_000;
        const double mean = 0.4;
        var arm = new TerminalArm(direction, armLen, armLen);
        var eventP = new CNEventPars(CNEventType.TelomereDuplication, 1, mean);
        var proportions = Enumerable.Range(0, 100000).Select(_ =>
        {
            var eventData = new TailEventData(_rnd, eventP, 0, contigLen, arm);
            long length = direction ? eventData.Start : contigLen - eventData.Start;
            return length / (double) armLen;
        });

        Assert.AreEqual(mean, proportions.Average(), 0.01);
    }

    [TestCase(true)]
    [TestCase(false)]
    public void TestTerminalEventLengthIsBoundedByTheArmProper(bool direction)
    {
        const long contigLen = 3_000_000;
        const long usableLen = 1_000_000;
        const long armLen = 1_100_000;
        // The Beta draw is taken against the fitted arm, which runs half a centromere further
        // than the arm proper, so every draw above the arm collapses onto a whole-arm event.
        var arm = new TerminalArm(direction, armLen, usableLen);
        var eventP = new CNEventPars(CNEventType.TelomereDeletion, 1, 0.6);
        int wholeArm = 0;

        for (int i = 0; i < 20000; i++)
        {
            var eventData = new TailEventData(_rnd, eventP, 0, contigLen, arm);
            long length = direction ? eventData.Start : contigLen - eventData.Start;
            Assert.LessOrEqual(length, usableLen);
            if (length == usableLen)
            {
                wholeArm++;
            }
        }

        // P(Beta(0.6, 0.4) > 1000000/1100000) is around a fifth of the draws.
        Assert.Greater(wholeArm, 0);
    }

    [Test]
    public void TestPairEventData()
    {
        const long lenA = 1_000_000;
        const long lenB = 10_000_000;
        var armA = new TerminalArm(true, 400_000, 400_000);
        var armB = new TerminalArm(false, 2_000_000, 2_000_000);
        var eventP = new CNEventPars(CNEventType.Translocation, 1, 0.1);
        var eventData = new PairEventData(
            _rnd, eventP, 0, lenA, armA, 1, lenB, armB);
        Assert.GreaterOrEqual(eventData.PosA, 0);
        Assert.Less(eventData.PosA, armA.UsableLength);
        Assert.Greater(eventData.PosB, 0);
        Assert.LessOrEqual(eventData.PosB, armB.UsableLength);
        Assert.IsTrue(eventData.Inverted);
    }

    [Test]
    public void TestBFBEventDataPlacesFinalBreakBetweenCentromeres()
    {
        const long len = 10_000_000;
        var arm = new TerminalArm(true, 4_000_000, 3_500_000);
        var eventP = new CNEventPars(CNEventType.BreakageFusionBridge, 1, 0.1);
        var eventData = new BFBEventData(_rnd, eventP, 0, len, arm);
        long removedLength = eventData.Start;
        long gapStart = len - arm.UsableLength;
        long gapEnd = len + arm.UsableLength - 2 * removedLength;
        Assert.GreaterOrEqual(removedLength, 1);
        Assert.LessOrEqual(removedLength, arm.UsableLength);
        Assert.GreaterOrEqual(eventData.FinalBreak, gapStart);
        Assert.LessOrEqual(eventData.FinalBreak, gapEnd);
    }

    [Test]
    public void TestCentromereBoundEventData()
    {
        const long len = 10_000_000;
        var cents = new List<(long start, long end)> { (1_000_000L, 2_000_000L) };
        var eventP = new CNEventPars(CNEventType.CentromereBoundDeletion, 1, 0.1);
        var eventData = new InternalEventData(_rnd, eventP, 0, cents, len);
        Assert.LessOrEqual(eventData.Start, cents[0].end);
        Assert.Less(eventData.Start, len);
        Assert.Greater(eventData.End, eventData.Start);
        Assert.Greater(eventData.End, cents[0].start);
    }
    
    [Test]
    public void TestArmEvent()
    {
        var eventP = new CNEventPars(CNEventType.ArmDeletion, 1);
        var cents = new List<(long start, long end)> { (1_000_000L, 2_000_000L) };
        var eventData = (TailEventData) new TailEventData(_rnd, eventP, 0, cents, 5_000_000L);
        Assert.Less(1_000_000L, eventData.Start);
        Assert.Greater(2_000_000L, eventData.Start);
    }

    [Test]
    public void TestTailEvent()
    {
        const long len = 10_000_000L;
        var eventP = new CNEventPars(CNEventType.TailDeletion, 1, 0.01);
        var eventData = new TailEventData(_rnd, eventP, 0, len);

        Assert.LessOrEqual(eventData.Length, len);
    }

    [TestCase(true)]
    [TestCase(false)]
    public void TestTelomereEventUsesPreselectedDirection(bool direction)
    {
        const long len = 10_000_000L;
        var eventP = new CNEventPars(CNEventType.TelomereDeletion, 1, 0.01);
        var eventData = new TailEventData(_rnd, eventP, 0, len, new TerminalArm(direction, len, len));

        Assert.AreEqual(direction, eventData.Direction);
        StringAssert.Contains(
            $"start:{(direction ? 0 : eventData.Start)};",
            eventData.EventDesc());
    }

    [Test]
    public void TestPyrgo()
    {
        const long len = 1_000_000;
        var eventP = new CNEventPars(CNEventType.Pyrgo, 1, 0.1, 10);
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
        var eventP = new CNEventPars(CNEventType.Rigma, 1, 0.1, 10);
        var eventData = new RigmaEventData(_rnd, eventP, 0, len);
        Assert.GreaterOrEqual(eventData.Start, 0);
        foreach (long stop in eventData.StopsList)
        {
            Assert.LessOrEqual(stop - eventData.Start, len);
        }
    }

    [Test]
    public void TestTemplatedEvent()
    {
        var eventP = new CNEventPars(CNEventType.TIBridge, 1, 0.1, 10);
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
    public void TestChromoplexy([Values] IntEdgeCases seed)
    {
        var eventP = new CNEventPars(CNEventType.Chromoplexy, 1, 0.1, 5);
        var sequences = new List<(int, long)>
        {
            (0, 5_000_000),
            (1, 10_000_000),
            (2, 20_000_000),
            (3, 30_000_000)
        };
        var eventData = new ChromoplexyEventData(new Random((int) seed), eventP, sequences);
        Assert.AreEqual(eventData.Stops.Sum(s => s.Count + 1), eventData.Sequence.Count);
        Console.WriteLine(eventData);
    }

    [Test]
    public void TestChromothripsis([Values] IntEdgeCases seed)
    {
        const long len = 100_000_000;
        var eventP = new CNEventPars(CNEventType.Chromothripsis, 1, 0.1);
        var eventData = new ChromothripsisEventData(new Random((int) seed), eventP, 0, len);
        Assert.Less(0, eventData.StopsList.Count);
        foreach (long stop in eventData.StopsList)
        {
            Assert.Greater(stop, 0);
            Assert.Less(stop, len);
        }
        Assert.GreaterOrEqual(eventData.StopsList.Count + 1,eventData.SelectionList.Count);
    }
    
    [Test]
    public void TestPointMutationEventData([Values] IntEdgeCases seed)
    {
        const long len = 100_000_000;
        var eventP = new CNEventPars(CNEventType.SNV, 1);
        var eventData = new PointMutationData(new Random((int) seed), eventP, 0, len);
        var dataArray = eventData.EventDesc().Split(';');
        Assert.AreEqual("contig:0", dataArray[0]);
        Assert.Less(eventData.Location, len);
        Assert.AreNotEqual("N", eventData.Base.ToString());
    }
}
