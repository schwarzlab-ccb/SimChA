using System.Collections;
using MathNet.Numerics.Distributions;
using SimChA.Computation;
using SimChA.DataTypes;

namespace SimChA.Simulation;

public class Simulator
{
    public List<List<SubClone>> Populations { get; }
    public IEnumerable<SubClone> FlatPops => CellSampling.Flatten(Populations);
    public SimParams SimParams { get; }

    private int newId;
    private int GetNewId() => ++newId;

    private int _generation;

    private double slowDownRate;

    private Random Rnd { get; }

    public Simulator(SimParams simParams, Random rnd)
    {
        SimParams = simParams;
        Rnd = rnd;
        var refKaryotype = new Karyotype(simParams.IsFemale, Rnd);
        var firstClone = new SubClone(0, -1, 0, SimParams.DivisionRate, refKaryotype, SimParams.InitialPop);
        Populations = new List<List<SubClone>> { new() { firstClone } };
    }

    public void Step()
    {
        _generation++;
        Kill();
        DivideAndMutate();
    }

    private void Kill()
    {
        foreach (var subClone in FlatPops.Where(sc => sc.AliveCount > 0))
        {
            int currentPop = subClone.AliveCount;
            int currentDead = subClone.DeadCount;
            int deadCount = Binomial.Sample(Rnd, SimParams.DeathRate, currentPop);
            subClone.NewGen(currentPop - deadCount, currentDead + deadCount);
        }
    }

    private bool IsViable(Karyotype kar) 
        => kar.ChromCount > 23;

    private void DivideAndMutate()
    {
        List<SubClone> newPops = new();
        foreach (var pop in Populations)
        {
            List<SubClone> newClones = new();
            foreach (var subClone in pop.Where(sc => sc.AliveCount > 0))
            {
                slowDownRate = Math.Pow(CellSampling.PopulationSize(pop), 1 / 3f) * SimParams.DivisionSlowDown;
                double divRate = subClone.DivisionRate - slowDownRate;
                if (divRate < 0)
                {
                    continue;
                }

                //  Cell division happens, split some of the cells
                int newCellsCount = Binomial.Sample(Rnd, divRate, subClone.AliveCount);
                int splitCellsCount = Binomial.Sample(Rnd, SimParams.SplitRate, newCellsCount);
                for (int splitI = 0; splitI < splitCellsCount; splitI++)
                {
                    var childClone = subClone.CreateChild(GetNewId(), _generation, 1.0);
                    newPops.Add(childClone);
                }

                // Mutate some of the cells
                int newMutantCount = Binomial.Sample(Rnd, SimParams.MutationRate, newCellsCount);
                int splitMutantCount = Binomial.Sample(Rnd, SimParams.SplitRate, newMutantCount);
                for (int mutationI = 0; mutationI < newMutantCount; mutationI++)
                {
                    double divChange = Rnd.NextDouble() < SimParams.DriverProb ? SimParams.FitnessInc : 1f;
                    var childClone = subClone.CreateChild(GetNewId(), _generation, divChange);
                    var abberation = SelectMutation();
                    childClone.Karyotype.ApplyAbberation(abberation);
                    if (!IsViable(childClone.Karyotype))
                    {
                        continue;
                    }

                    // Split some of the mutatnts
                    if (mutationI < splitMutantCount)
                    {
                        newPops.Add(childClone);
                    }
                    else
                    {
                        newClones.Add(childClone);
                    }
                }

                subClone.AddNewCells(newCellsCount - splitCellsCount - newMutantCount);
            }

            pop.AddRange(newClones);
        }

        // Create new population from the split cells
        Populations.AddRange(newPops.Select(sc => new List<SubClone> { sc }));
    }

    private AbberationEnum SelectMutation()
    {
        double ratesSum = SimParams.RatesSum;
        double sample = ContinuousUniform.Sample(Rnd, 0, ratesSum);
        foreach ((var abb, double rate) in SimParams.AbberationRates)
        {
            if (sample <= rate)
            {
                return abb;
            }
            sample -= rate;
        }

        // In case float-point calculations would cause jumping out of the loop
        return SimParams.AbberationRates.Last().Key;
    }
}