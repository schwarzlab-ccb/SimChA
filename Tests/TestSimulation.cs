// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using System;
using MathNet.Numerics.Distributions;
using NUnit.Framework;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.Simulation;

namespace Tests;

[TestFixture]
public class TestSimulation
{
    private SimParams _simParams;
    private Simulator _sim;
    
    [SetUp]
    public void Setup()
    {
        _simParams = new SimParams
        {
            DivisionRate = 0.1f, 
            MutationRate = 0.01f, 
            IsFemale = true,
            AbberationRates =
            {
                [AbberationEnum.InternalDeletion] = 50f,
                [AbberationEnum.InternalDuplication] = 50f,
                [AbberationEnum.Translocation] = 20f,
                [AbberationEnum.TailDeletion] = 15f,
                [AbberationEnum.BreakageFusionBridge] = 10f,
                [AbberationEnum.Inversion] = 10f,
                [AbberationEnum.Missegregation] = 5f,
                [AbberationEnum.Duplication] = 5f,
                [AbberationEnum.Chromothripsis] = 1f
            }
        };
        _sim = new Simulator(_simParams, new Random(0));
    }

    [Test]
    public void TestBinomial()
    {
        var rnd1 = new Random(0);
        var rnd2 = new Random(0);
        int res1 = 0, res2 = 0;
        for (int i = 0; i < 100; i++)
        {
            res1 = Binomial.Sample(rnd1, .5f, 1000);
        }
        for (int i = 0; i < 100; i++)
        {
            res2 = Binomial.Sample(rnd2, .5f, 1000);
        }
        Assert.AreEqual(res1, res2);
    }
}