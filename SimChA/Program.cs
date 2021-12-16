using MathNet.Numerics.Distributions;
using SimChA.DataTypes;
using SimChA.Simulation;

Console.WriteLine("SimChA 0.0.1");

var simParams = new SimParams { DivisionRate = .1f, MutationRate = 0.01f, IsFemale = true};
var simulator = new Simulator(simParams);
for (int i = 0; i < 150; i++)
{
    simulator.Step();
}
Console.WriteLine($"Cell count {simulator.Clones.Sum(c => c.AliveCount)}");
Console.WriteLine($"SubClone count {simulator.Clones.Count}");
Console.WriteLine($"Last Karyotype {simulator.Clones.Last().Karyotype}");