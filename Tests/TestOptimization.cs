// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SimChA.Optimization;

namespace Tests;

public class TestOptimization
{
    [Test]
    public void TestWD()
    {
        var a = new List<int>() { 0, 1, 3 };
        var b = new List<int>() { 1, 2, 1 };
        var aCDF = StatisticMeasures<int>.GetCDF(a);
        Assert.AreEqual(.25, aCDF[1], double.Epsilon);
        Assert.AreEqual(1, aCDF[2], double.Epsilon);
        var dist = StatisticMeasures<int>.WassersteinDistance(a, b);
    }
}