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
    public void TestWD2()
    {
        var a = new List<double>() {0, 4, 0, 0};
        var b = new List<double>() {2, 0, 1, 1};
        var aCDF = StatisticMeasures.GetCDF(a);
        var bCDF = StatisticMeasures.GetCDF(b);
        var dist = StatisticMeasures.WassersteinDistance(a, b);
        Assert.AreEqual(1, aCDF[2], double.Epsilon);
    }

    [Test]
    public void TestWDHistogram()
    {
        var a = new List<double>() {1, 1, 1, 1};
        var b = new List<double>() {0, 2, 3, 0};
        var max = Math.Max(a.Max(), b.Max());
        var min = Math.Min(a.Min(), b.Min());
        var bins = 4;
        var aHist = new Histogram(a, bins, min, max);
        var bHist = new Histogram(b, bins, min, max);
        Assert.AreEqual(0, aHist[0].Count, double.Epsilon);
        Assert.AreEqual(4, aHist[1].Count, double.Epsilon);
        Assert.AreEqual(0, aHist[2].Count, double.Epsilon);
        Assert.AreEqual(0, aHist[3].Count, double.Epsilon);

        Assert.AreEqual(4, aHist.DataCount, double.Epsilon);

        Assert.AreEqual(2, bHist[0].Count, double.Epsilon);
        Assert.AreEqual(0, bHist[1].Count, double.Epsilon);
        Assert.AreEqual(1, bHist[2].Count, double.Epsilon);
        Assert.AreEqual(1, bHist[3].Count, double.Epsilon);
        Assert.AreEqual(4, bHist.DataCount, double.Epsilon);
        var aCDF = StatisticMeasures.GetCDF(aHist);
        Assert.AreEqual(1, aCDF[1], double.Epsilon);
        var bCDF = StatisticMeasures.GetCDF(bHist);
        Assert.AreEqual(0.5, bCDF[0], double.Epsilon);
        var dist = StatisticMeasures.WassersteinDistance(aHist, bHist);
        var binWidth = (max - min )/ bins;
        Assert.AreEqual(0.3125, dist, double.Epsilon);
    }
}