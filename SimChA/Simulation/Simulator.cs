using MathNet.Numerics.Distributions;
using SimChA.Computation;
using SimChA.DataTypes;

namespace SimChA.Simulation;

public class Simulator
{
    public List<List<SubClone>> Populations { get; }
    public IEnumerable<SubClone> FlatPops => CellSampling.Flatten(Populations);
    public SimParams SimParams { get; }

    public ProgramConfig Config { get; }

    public int AliveSC;

    private int newId;
    private int GetNewId() => ++newId;

    private int _generation;
    private Random Rnd { get; }

    public Simulator(SimParams simParams, ProgramConfig config, Random rnd)
    {
        SimParams = simParams;
        Config = config;
        Rnd = rnd;
        // var refKaryotype = new Karyotype(simParams.IsFemale, Rnd);
        // int popSize = (int)Math.Round(1 / SimParams.MutationRate);
        // popSize = simParams.InitialPop;
        var firstClone = new SubClone(0, -1, 0, SimParams.DivisionRate * (1 + SimParams.FitnessLambdaInv), 1, simParams.InitialPop);
        Populations = new List<List<SubClone>> { new() { firstClone } };
    }

    public void Step()
    {
        _generation++;
        List<SubClone> newPops = new();
        foreach (var pop in Populations)
        {
            List<SubClone> newClones = new();
            long popSize = CellSampling.PopulationSize(pop);
            long aliveCount = CellSampling.AliveCount(pop);
            double divisionFraction = 1f;
            if (SimParams.Confinement > 0 && aliveCount > 1 / SimParams.Confinement)
            {
                double divisible = Math.Round(Math.Pow(popSize, 2 / 3f)) / SimParams.Confinement;
                if (aliveCount > divisible && aliveCount > 0)
                {
                    divisionFraction = divisible / aliveCount;
                }
            }

            AliveSC = 0;
            foreach (var subClone in pop.Where(sc => sc.AliveCount > 0))
            {
                AliveSC++;

                // Kill cells
                int newDead;
                if (Config.StochasticCellLife)
                {
                    newDead = Binomial.Sample(Rnd, SimParams.DivisionRate, (int) subClone.AliveCount);
                }
                else
                {
                    double deadFraction = SimParams.DivisionRate * subClone.AliveCount + subClone.ToDie;
                    newDead = (int)deadFraction;
                    subClone.ToDie = deadFraction - newDead;
                }
                
                // Create new cells
                int newCellsCount;
                double divRate = Math.Clamp(subClone.DivisionRate * divisionFraction, 0.0, 1.0);
                if (Config.StochasticCellLife)
                {
                    newCellsCount = Binomial.Sample(Rnd, divRate, (int) subClone.AliveCount);
                }
                else
                {
                    double divideFraction = (divRate > 0 ? divRate * subClone.AliveCount : 0) + subClone.ToDivide;
                    newCellsCount = (int)divideFraction;
                    subClone.ToDivide = divideFraction - newCellsCount;
                }

                // Decayed cells
                // int newDecayed = SimParams.DecayRate > 0 ? Binomial.Sample(Rnd, SimParams.DecayRate, subClone.DeadCount) : 0;


                //  From some of the cells, create new populations
                int splitCellsCount = Binomial.Sample(Rnd, SimParams.SplitRate, newCellsCount);
                for (int splitI = 0; splitI < splitCellsCount; splitI++)
                {
                    var childClone = subClone.CreateChild(
                        GetNewId(), _generation, subClone.DivisionRate, subClone.NumberDrivers);
                    newPops.Add(childClone);
                }

                // Mutate some of the cells
                int newMutantCount = SimParams.MutationRate > 0 
                    ? Binomial.Sample(Rnd, SimParams.MutationRate, newCellsCount * 2) 
                    : 0;
                int splitMutantCount = SimParams.SplitRate > 0
                    ? Binomial.Sample(Rnd, SimParams.SplitRate, newMutantCount)
                    : 0;
                for (int mutationI = 0; mutationI < newMutantCount; mutationI++)
                {
                    double divChange = Math.Min(
                        Exponential.Sample(Rnd, 1/SimParams.FitnessLambdaInv) * SimParams.DivisionRate,
                        3 * SimParams.FitnessLambdaInv);

                    double newDivision = Config.MultiplicativeFitness
                        ? subClone.DivisionRate * divChange
                        : subClone.DivisionRate + divChange;
                    var childClone = subClone.CreateChild(GetNewId(), _generation, newDivision, subClone.NumberDrivers + 1);

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
                    (uint) (subClone.AliveCount + newCellsCount - splitCellsCount - newMutantCount - newDead),
                    (uint) (subClone.DeadCount + newDead));
            }

            pop.AddRange(newClones);
        }

        // Create new population from the split cells
        Populations.AddRange(newPops.Select(sc => new List<SubClone> { sc }));
    }

    // private bool IsViable(Karyotype kar)
    //     => kar.ChromCount > 23;
    //
    // private AberrationEnum SelectMutation()
    // {
    //     double ratesSum = SimParams.SumRates();
    //     double sample = ContinuousUniform.Sample(Rnd, 0, ratesSum);
    //     foreach ((var abb, double rate) in SimParams.AberrationRates)
    //     {
    //         if (sample <= rate)
    //         {
    //             return abb;
    //         }
    //
    //         sample -= rate;
    //     }
    //
    //     // In case float-point calculations would cause jumping out of the loop
    //     return SimParams.AberrationRates.Last().Key;
    // }
}