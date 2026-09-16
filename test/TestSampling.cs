using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Statistics;
using NUnit.Framework;
using SimChA.Computation;
using SimChA.Simulation;

namespace Tests;

[TestFixture]
public class TestSampling
{
    private Random _rnd = null!;
    
    [SetUp]
    public void Setup()
    {
        _rnd = new Random(0);
    }

    // The interior length distribution. The predecessor asserted only that the realized mean landed
    // within 50% of the configured one, which is why a systematic 10-13% shortfall went unnoticed
    // for as long as it did: GetParetoScale had no shape term and so inverted no mean at all.
    [Test]
    public void TestParetoRecoversMean(
        [Values(0.005, 0.05, 0.1058, 0.1312, 0.3)] double mean,
        [Values(0.12, 0.378, 0.5, 1.0, 2.5)] double shape)
    {
        double sampled = Enumerable.Range(0, 200_000)
            .Select(_ => Sampling.SampleParetoLim(_rnd, mean, shape)).Mean();
        Assert.AreEqual(mean, sampled, 0.004 + mean * 0.03);
    }

    [Test]
    public void TestParetoShapeDefaultsToTheFixedValue()
    {
        Assert.AreEqual(0.5, Sampling.FixedParetoShape);
        var implicitShape = Enumerable.Range(0, 500)
            .Select(_ => Sampling.SampleParetoLim(new Random(11), 0.1058)).ToList();
        var explicitShape = Enumerable.Range(0, 500)
            .Select(_ => Sampling.SampleParetoLim(new Random(11), 0.1058, Sampling.FixedParetoShape))
            .ToList();
        Assert.AreEqual(implicitShape, explicitShape);
    }

    // The scale is what absorbs the shape: holding the mean fixed, a heavier tail has to start
    // lower. If this ever stops depending on the shape, the old shape-free formula is back.
    [Test]
    public void TestParetoScaleSolvesTheConfiguredMean(
        [Values(0.05, 0.1312, 0.3)] double mean, [Values(0.12, 0.5, 2.5)] double shape)
    {
        double scale = Sampling.GetParetoScale(mean, shape);
        Assert.AreEqual(mean, Sampling.BoundedParetoMean(scale, shape), 1e-9);
    }

    [Test]
    public void TestLowerParetoShapeIsHeavierTailed()
    {
        // Both draws have the same mean by construction, so the shape shows up only in the spread:
        // a heavier tail carries the mean further out and leaves the body shorter. This is the
        // degree of freedom the empirical fit wanted -- interior gain 0.378, interior loss 0.121.
        var heavy = Enumerable.Range(0, 200_000)
            .Select(_ => Sampling.SampleParetoLim(_rnd, 0.1312, 0.12)).ToList();
        var light = Enumerable.Range(0, 200_000)
            .Select(_ => Sampling.SampleParetoLim(_rnd, 0.1312, 0.5)).ToList();
        Assert.AreEqual(heavy.Mean(), light.Mean(), 0.01);
        Assert.Greater(heavy.Quantile(0.9), light.Quantile(0.9));
        Assert.Less(heavy.Median(), light.Median());
    }

    [Test]
    public void TestParetoSupportIsTheUnitInterval(
        [Values(0.01, 0.1312, 0.4)] double mean, [Values(0.05, 0.5, 4.0)] double shape)
    {
        var samples = Enumerable.Range(0, 50_000)
            .Select(_ => Sampling.SampleParetoLim(_rnd, mean, shape)).ToList();
        Assert.GreaterOrEqual(samples.Min(), 0.0);
        Assert.LessOrEqual(samples.Max(), 1.0);
    }

    [Test]
    public void TestParetoRejectsUnreachableParameters()
    {
        // The bounded distribution tends to uniform on [0,1] as its scale grows, so one half is a
        // supremum no shape can reach; asking for it is a configuration error, not something to
        // silently approximate.
        Assert.Throws<ArgumentOutOfRangeException>(() => Sampling.GetParetoScale(0.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => Sampling.GetParetoScale(0.8));
        Assert.Throws<ArgumentOutOfRangeException>(() => Sampling.GetParetoScale(0.1, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => Sampling.GetParetoScale(0.1, -1));
    }

    [Test]
    public void TestBoundedParetoMeanIsNumericallyStable()
    {
        // The bounded distribution flattens to uniform on [0,1] as the scale grows, so the mean
        // approaches one half from below and never passes it. The direct form of the integral
        // cancels to noise past about 1e6 -- it returned 828 at L = 1e9 and 1.9e9 at L = 1e12.
        foreach (double shape in new[] { 0.1, 0.5, 1.0, 3.0 })
        {
            // Well inside the range the solver searches, where the value must be exact.
            Assert.AreEqual(0.4999875, Sampling.BoundedParetoMean(1e4, 0.5), 1e-6);
            // Far outside it, where only the absence of a blow-up is being claimed.
            foreach (double scale in new[] { 1e6, 1e9, 1e12 })
            {
                Assert.AreEqual(0.5, Sampling.BoundedParetoMean(scale, shape), 2e-3);
            }
            // The smallest scale searched drives log((L+1)/L) to 691, which overflows the
            // expm1 form; below one the direct form is used instead and stays finite.
            double tiny = Sampling.BoundedParetoMean(1e-300, shape);
            Assert.That(tiny, Is.GreaterThanOrEqualTo(0).And.LessThan(1e-3));
        }
    }

    [Test]
    public void TestParetoSegNeedsNoConfiguredMean()
    {
        // Event types that carry no Frac reach this path in tests; one base is the same floor
        // GetBetaSeg applies, rather than a solve for an unreachable mean.
        Assert.AreEqual(1, Sampling.GetParetoSeg(_rnd, 1_000_000, 0));
    }

    // The telomere-bound length distribution. Unlike the Pareto side this was always exactly
    // parameterised -- beta = alpha*(1-mean)/mean gives Beta(alpha, beta) the requested mean for
    // any alpha -- so configuring alpha changes the spread and leaves the mean alone.
    [Test]
    public void TestBetaSamplingRecoversMean(
        [Values(0.25, 0.5, 0.5954, 0.75)] double mean, [Values(0.3, 0.6, 1.5)] double alpha)
    {
        double sampledMean = Enumerable.Range(0, 100000)
            .Select(_ => Sampling.SampleBeta(_rnd, mean, alpha))
            .Average();
        Assert.AreEqual(mean, sampledMean, 0.01);
    }

    [Test]
    public void TestBetaAlphaDefaultsToTheFixedValue()
    {
        var implicitAlpha = Enumerable.Range(0, 500)
            .Select(_ => Sampling.SampleBeta(new Random(13), 0.5954)).ToList();
        var explicitAlpha = Enumerable.Range(0, 500)
            .Select(_ => Sampling.SampleBeta(new Random(13), 0.5954, Sampling.FixedBetaAlpha))
            .ToList();
        Assert.AreEqual(implicitAlpha, explicitAlpha);
    }

    [Test]
    public void TestLowerBetaAlphaIsMoreUShaped()
    {
        // Smaller alpha piles mass at both ends without moving the mean, which is what a fit of the
        // telomere-bound class would adjust.
        var spread = Enumerable.Range(0, 200_000)
            .Select(_ => Sampling.SampleBeta(_rnd, 0.5954, 0.25)).ToList();
        var tight = Enumerable.Range(0, 200_000)
            .Select(_ => Sampling.SampleBeta(_rnd, 0.5954, 1.5)).ToList();
        Assert.AreEqual(spread.Mean(), tight.Mean(), 0.01);
        Assert.Greater(spread.StandardDeviation(), tight.StandardDeviation());
    }

    [Test]
    public void TestBetaShapeIsDerivedFromMean()
    {
        Assert.AreEqual(0.6, Sampling.FixedBetaAlpha);
        Assert.AreEqual(0.6, Sampling.GetBetaShape(0.5), 1e-12);
        Assert.AreEqual(1.8, Sampling.GetBetaShape(0.25), 1e-12);
        Assert.AreEqual(0.25, Sampling.GetBetaShape(0.5, 0.25), 1e-12);
        Assert.Throws<ArgumentOutOfRangeException>(() => Sampling.GetBetaShape(1));
        Assert.Throws<ArgumentOutOfRangeException>(() => Sampling.GetBetaShape(0.5, 0));
    }

    [Test]
    public void TestContSampling([Values(DistType.Exponential, DistType.Normal)] DistType dist, [Values(0.01, 0.1, 1, 10, 100)] double mean)
    {
        double res = Enumerable.Range(0, 1000).Select(i => Sampling.SampleContDist(_rnd, dist, mean)).Average();
        Assert.Greater(res, mean * 0.5);
        Assert.Less(res, mean * 1.5);
    }
    
    [Test]
    public void TestDiscSampling([Values(DistType.Geometric, DistType.Poisson)] DistType dist, [Values(1, 10, 100)] double mean)
    {
        double res = Enumerable.Range(0, 1000).Select(i => Sampling.SampleDiscDist(_rnd, dist, mean)).Average();
        Assert.Greater(res, mean * 0.5);
        Assert.Less(res, mean * 1.5);
    }

    [Test]
    public void TestGeometric()
    {
        double res = Enumerable.Range(0, 1000000).Select(i => Sampling.SampleDiscDist(_rnd, DistType.Geometric, 62)).Average();
        Assert.AreEqual(res, 62, 0.1);
    }

    [Test]
    public void TestGammaRecoversMeanAndShape([Values(1.5, 4.0, 13.6)] double shape)
    {
        var samples = Enumerable.Range(0, 1000000)
            .Select(i => (double) Sampling.SampleDiscDist(_rnd, DistType.Gamma, 62, shape)).ToList();
        // Shape sets the dispersion only: mean must come out at RateMean regardless of it, and
        // the CV at 1/sqrt(shape) -- the property that lets one shape serve many cancer types.
        Assert.AreEqual(62, samples.Mean(), 0.1);
        Assert.AreEqual(1 / Math.Sqrt(shape), samples.StandardDeviation() / samples.Mean(), 0.01);
    }

    [Test]
    public void TestGammaAtShapeOneMatchesGeometric()
    {
        // shape 1 is the exponential, so it should reproduce the geometric SimChA sampled before.
        var gamma = Enumerable.Range(0, 200000)
            .Select(i => (double) Sampling.SampleDiscDist(_rnd, DistType.Gamma, 62, 1)).ToList();
        Assert.AreEqual(62, gamma.Mean(), 0.5);
        Assert.AreEqual(1.0, gamma.StandardDeviation() / gamma.Mean(), 0.02);
    }

    [Test]
    public void TestGammaCountsAreUsable()
    {
        // A tiny mean must still yield at least one event rather than a zero-event sample.
        var counts = Enumerable.Range(0, 10000)
            .Select(i => Sampling.SampleDiscDist(_rnd, DistType.Gamma, 0.5, 1.5)).ToList();
        Assert.GreaterOrEqual(counts.Min(), 1);
    }

    [Test]
    public void TestGammaIsSeedDeterministic()
    {
        var first = Enumerable.Range(0, 100)
            .Select(i => Sampling.SampleDiscDist(new Random(7), DistType.Gamma, 62, 1.5)).ToList();
        var second = Enumerable.Range(0, 100)
            .Select(i => Sampling.SampleDiscDist(new Random(7), DistType.Gamma, 62, 1.5)).ToList();
        Assert.AreEqual(first, second);
    }

    [Test]
    public void TestStops()
    {
        var stops = Sampling.GetStopsForShards(_rnd, 1_000_000, 10);
        Assert.AreEqual(9, stops.Count);
        long prev = 0;
        foreach (long stop in stops)
        {
            Assert.Greater(stop, 0);
            Assert.Less(stop, 1_000_000);
            Assert.Greater(stop, prev);
            prev = stop;
        }
    }
    
    [Test]
    public void TestGetExpSeg()
    {
        Assert.AreEqual(31987, Sampling.GetExpSeg(_rnd, 1_000_000, 0.1));
        Assert.AreEqual(202, Sampling.GetExpSeg(_rnd, 1_000_000, 1_000));
    }
    
    [Test]
    public void TestGetNormSeg()
    {
        Assert.AreEqual(119319, Sampling.GetNormSeg(_rnd, 1_000_000, 0.1));
    }

    [Test]
    public void TestDirichlet()
    {
        var mixture = Sampling.CreateRandomMixture(_rnd, [1000.0, 1.0, 0.001]);
        Assert.AreEqual(3, mixture.Count);
        Assert.AreEqual(1.0, mixture.Sum(), 1e-6);
        mixture.ForEach(x => Assert.GreaterOrEqual(x, 0.0));
        Assert.GreaterOrEqual(mixture[0], mixture[1]);
        Assert.GreaterOrEqual(mixture[1], mixture[2]);
    }
    
    [Test]
    public void TestIndexPicker()
    {
        for (int i = 0; i < 10; i++)
        {
            var rnd = new Random();
            var probs = new List<double> { 0.1, 0.2, 0.0, -1.0, 0.3, 0.4 };
            int index = rnd.PickRndIndex(probs);
            Assert.GreaterOrEqual(index, 0);
            Assert.Less(index, probs.Count);
            Assert.AreNotEqual(index, 2);
            Assert.AreNotEqual(index, 3);
        }
    }
}
