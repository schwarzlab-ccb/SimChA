// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using System;
using System.Text.Json;
using NUnit.Framework;
using SimChA.IO;
using SimChA.Simulation;

namespace Tests;

[TestFixture]
public class TestIO
{
    [Test]
    public void TestConfigSerialization()
    {
        var defaultAberrations = Aberrations.DefaultAberrations();
        var simParams = SimParams.CreateSimParams(0, true, defaultAberrations);
        var options = new JsonSerializerOptions { WriteIndented = true };
        string serialized = JsonSerializer.Serialize(simParams, options);
        Console.WriteLine(serialized);
        var deserialized = JsonSerializer.Deserialize<SimParams>(serialized);
        Assert.NotNull(deserialized);
        Assert.AreEqual(simParams, deserialized);
    }
}