using SimChA.Computation;
using SimChA.DataTypes;
using static Extreme.Statistics.Distributions.BinomialDistribution;
using ExtremeBinDist = Extreme.Statistics.Distributions.BinomialDistribution;

namespace SimChA.Simulation;

public class Simulator
{
    public int AliveClones => Clones.Count(c => c.IsAlive);
    public int StepNo;
    private int LastId;

    public Simulator(SimParams simParams, Random rnd)
    {
        Rnd = rnd;
        StepNo = 0;
        LastId = -1;
        SimParams = simParams;
        var primeval = new SubClone(GetNewId(), -1, SimParams.StartMut, SimParams.StartPop);
        Clones = new List<SubClone> { primeval };
    }

    public List<SubClone> Clones { get; }
    public long CellCount => Clones.Sum(clone => clone.CellCount);
    public SimParams SimParams { get; }
    private Random Rnd { get; }

    private int GetNewId() => ++LastId;

    public void Step()
    {
        StepNo++;

        List<SubClone> newClones = new();

        foreach (var clone in Clones.Where(c => c.IsAlive))
        {
            int newDead = Sample(Rnd, clone.CellCount, SimParams.DeathRate * SimParams.Turnover);
            int divisions = Sample(Rnd, clone.CellCount, SimParams.Turnover);
            int mutations = Sample(Rnd, 2 * divisions, SimParams.MutationProb);
            clone.CellCount = Math.Max(0, clone.CellCount + divisions - newDead - mutations);
            for (int i = 0; i < mutations; i++)
            {
                newClones.Add(clone.CreateChild(GetNewId()));
            }
        }

        Clones.AddRange(newClones);
    }
}