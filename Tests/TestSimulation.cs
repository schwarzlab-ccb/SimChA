// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using System;
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
        _sim = new Simulator(_simParams);
    }

    [Test]
    public void TestSampling()
    {
        while (_sim.Clones.Count < 100)
        {
            _sim.Step();
        }
        var sampleCells = CellSampling.SampleCells(_sim.Clones, 1000);
        Console.WriteLine($"Pop size: {CellSampling.PopulationSize(sampleCells)}");
    }
}