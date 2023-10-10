// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Statistics;
using NUnit.Framework;
using SimChA.Optimization;

namespace Tests;

public class TestOptimization
{
    [Test]
    public void TestWD()
    {
        var a = new List<double>() { 0, 1, 3 };
        var b = new List<double>() { 1, 2, 1 };
        var aCDF = StatisticMeasures.GetCDF(a);
        Assert.AreEqual(.25, aCDF[1], double.Epsilon);
        Assert.AreEqual(1, aCDF[2], double.Epsilon);
        var dist = StatisticMeasures.WassersteinDistance(a, b);
        Assert.AreEqual(0.25, dist, double.Epsilon);
    }

    [Test]
    public void TestWDHistogram()
    {
        var a = new List<double>() {1, 1, 2};
        var b = new List<double>() {0, 3, 1};
        var max = Math.Max(a.Max(), b.Max());
        var min = Math.Min(a.Min(), b.Min());
        var bins = 4;
        var aHist = new Histogram(a, bins, min, max);
        var bHist = new Histogram(b, bins, min, max);
        Assert.AreEqual(0, aHist[0].Count, double.Epsilon);
        Assert.AreEqual(2, aHist[1].Count, double.Epsilon);
        Assert.AreEqual(1, aHist[2].Count, double.Epsilon);
        Assert.AreEqual(0, aHist[3].Count, double.Epsilon);

        Assert.AreEqual(3, aHist.DataCount, double.Epsilon);

        Assert.AreEqual(1, bHist[0].Count, double.Epsilon);
        Assert.AreEqual(1, bHist[1].Count, double.Epsilon);
        Assert.AreEqual(0, bHist[2].Count, double.Epsilon);
        Assert.AreEqual(1, bHist[3].Count, double.Epsilon);
        var aCDF = StatisticMeasures.GetCDF(aHist);
        Assert.AreEqual(2.0/3, aCDF[1], double.Epsilon);
        var bCDF = StatisticMeasures.GetCDF(bHist);
        Assert.AreEqual(1.0/3, bCDF[0], double.Epsilon);
        var dist = StatisticMeasures.WassersteinDistance(aHist, bHist);
        Assert.AreEqual(1.0/6, dist, double.Epsilon);
    }
}