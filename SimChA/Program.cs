using SimChA.DataTypes;
using SimChA.Simulation;

Console.WriteLine("SimChA 0.0.1");

var simParams = new SimParams
{
    DivisionRate = .1f, MutationRate = 0.01f, IsFemale = true,
    AbberationRates =
    {
        [AbberationEnum.InternalDeletion] = 100f,
        [AbberationEnum.Duplication] = 50f,
        [AbberationEnum.Translocation] = 20f,
        [AbberationEnum.TailDeletion] = 15f,
        [AbberationEnum.BreakageFusionBridge] = 10f,
        [AbberationEnum.Inversion] = 10f,
        [AbberationEnum.Missegregation] = 5f,
        [AbberationEnum.Chromothripsis] = 1f
    }
};

var simulator = new Simulator(simParams);
const int simSteps = 150;
for (int i = 0; i < simSteps; i++)
{
    simulator.Step();
    Console.WriteLine($"Steps {i + 1}/{simSteps}");
}

Console.WriteLine($"Cell count {simulator.Clones.Sum(c => c.AliveCount)}");
Console.WriteLine($"SubClone count {simulator.Clones.Count}");
Console.WriteLine($"Last Karyotype: {simulator.Clones.Last().Karyotype}");