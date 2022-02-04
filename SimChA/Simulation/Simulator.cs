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
        List<SubClone> newPops = new();
        foreach (var pop in Populations)
        {
            List<SubClone> newClones = new();
            slowDownRate = Math.Pow(CellSampling.AliveCount(pop), 1f / 3) * SimParams.DivisionSlowDown;
            
            foreach (var subClone in pop.Where(sc => sc.AliveCount > 0))
            {
                // Kill cells
                int newDead = Binomial.Sample(Rnd, SimParams.DeathRate, subClone.AliveCount);
                
                // Create new cells
                double divRate = Math.Clamp(subClone.DivisionRate - slowDownRate, 0, 1);
                int newCellsCount = Binomial.Sample(Rnd, divRate, subClone.AliveCount);

                //  From some of the cells, create new populations
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
                    var aberration = SelectMutation();
                    childClone.Karyotype.ApplyAbberation(aberration);
                    if (!IsViable(childClone.Karyotype))
                    {
                        continue;
                    }
                    // Some of the mutations may create new populations
                    if (mutationI < splitMutantCount)
                    {
                        newPops.Add(childClone);
                    }
                    else
                    {
                        newClones.Add(childClone);
                    }
                }
                subClone.NewGen(
                    subClone.AliveCount + newCellsCount - splitCellsCount - newMutantCount - newDead,
                    subClone.DeadCount + newDead);
            }
            pop.AddRange(newClones);
        }

        // Create new population from the split cells
        Populations.AddRange(newPops.Select(sc => new List<SubClone> { sc }));
    }

    private bool IsViable(Karyotype kar) 
        => kar.ChromCount > 23;

    private AberrationEnum SelectMutation()
    {
        double ratesSum = SimParams.SumRates();
        double sample = ContinuousUniform.Sample(Rnd, 0, ratesSum);
        foreach ((var abb, double rate) in SimParams.AberrationRates)
        {
            if (sample <= rate)
            {
                return abb;
            }
            sample -= rate;
        }

        // In case float-point calculations would cause jumping out of the loop
        return SimParams.AberrationRates.Last().Key;
    }
}