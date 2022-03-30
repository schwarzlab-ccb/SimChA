using MathNet.Numerics.Distributions;
using SimChA.Computation;
using SimChA.DataTypes;
using  ExtremeBinDist = Extreme.Statistics.Distributions.BinomialDistribution;

namespace SimChA.Simulation;

public class Simulator
{
    public List<List<SubClone>> Populations { get; }
    public IEnumerable<SubClone> FlatPops => CellSampling.Flatten(Populations);
    public SimParams SimParams { get; }

    public int AliveSC;

    private int newId;
    private int GetNewId() => ++newId;

    private int _generation;
    private Random Rnd { get; }

    private const int INIT_POP = 1;
    const double THIRD = 1.0 / 3.0;
    const double TWO_THIRD = 2.0 / 3.0;

    public Simulator(SimParams simParams, Random rnd)
    {
        SimParams = simParams;
        Rnd = rnd;
        // var refKaryotype = new Karyotype(simParams.IsFemale, Rnd);
        // int popSize = (int)Math.Round(1 / SimParams.MutationRate);
        // popSize = simParams.InitialPop;
        var firstClone = new SubClone(0, -1, 0, BumpFitness(SimParams.BirthRate, SimParams, Rnd), 1, SimParams.InitPop);
        Populations = new List<List<SubClone>> { new() { firstClone } };
    }

    private static double BumpFitness(double original, SimParams simParams, Random rnd)
    {
        double divChange = FitnessFunction.SampleFitness(simParams, rnd);
        return simParams.MultiplicativeFitness ? original * divChange : original + divChange;
    }
    
    public void Step()
    {
        _generation++;
        AliveSC = 0;
        List<SubClone> newPops = new();
        foreach (var pop in Populations)
        {
            List<SubClone> newClones = new();
            long popSize = CellSampling.PopulationSize(pop);
            long aliveCount = CellSampling.AliveCount(pop);
            double divisionFraction = 1;
            if (SimParams.Confinement > 0 && aliveCount > 1 / SimParams.Confinement)
            {
                double divisible = Math.Round(Math.Pow(popSize, TWO_THIRD)) / SimParams.Confinement;
                if (aliveCount > divisible && aliveCount > 0)
                {
                    divisionFraction = Math.Clamp(divisible / aliveCount, 0.0, 1.0);
                }
            }

            foreach (var subClone in pop.Where(sc => sc.AliveCount > 0))
            {
                AliveSC++;

                // Kill cells
                int newDead;
                if (SimParams.StochasticCellLife)
                {
                    // newDead = Binomial.Sample(Rnd, SimParams.DivisionRate, (int)subClone.AliveCount);
                    newDead = ExtremeBinDist.Sample(Rnd, (int)subClone.AliveCount, SimParams.BirthRate);
                }
                else
                {
                    double deadFraction = SimParams.BirthRate * subClone.AliveCount + subClone.ToDie;
                    newDead = (int)deadFraction;
                    subClone.ToDie = deadFraction - newDead;
                }

                // Create new cells
                int newCellsCount;
                double divRate = Math.Clamp(subClone.DivisionRate * divisionFraction, 0.0, 1.0);
                if (SimParams.StochasticCellLife)
                {
                    newCellsCount =
                        ExtremeBinDist.Sample(Rnd, (int)subClone.AliveCount, divRate);
                    // newCellsCount = Binomial.Sample(Rnd, divRate, (int)subClone.AliveCount);
                }
                else
                {
                    double divideFraction = (divRate > 0 ? divRate * subClone.AliveCount : 0) + subClone.ToDivide;
                    newCellsCount = (int)divideFraction;
                    subClone.ToDivide = divideFraction - newCellsCount;
                }
                
                // Calculate Split chance 
                double splitFactor = SimParams.SplitRate > 0 ?
                    Math.Clamp((SimParams.SplitRate * Math.Pow(popSize, THIRD * (1f - divisionFraction)) / Populations.Count), 0.0, 1.0) : 0.0;

                // Mutate some of the cells
                int newMutantCount = SimParams.MutationRate > 0 
                    ?  ExtremeBinDist.Sample(Rnd, newCellsCount, SimParams.MutationRate) : 0;

                for (int mutationI = 0; mutationI < newMutantCount; mutationI++)
                {

                    double newDivision = BumpFitness(subClone.DivisionRate, SimParams, Rnd);
                    var childClone = subClone.CreateChild(GetNewId(), _generation, newDivision, subClone.NumberDrivers + 1);
                    if (splitFactor > 0 && Rnd.NextDouble() < splitFactor)
                    {
                        newPops.Add(childClone);
                    }
                    else
                    {
                        newClones.Add(childClone);
                    }
                }
                
                //  From some of the cells, create new populations
                int splitCellsCount = 0;
                if (splitFactor > 0 && newCellsCount > newMutantCount)
                {
                    splitCellsCount = ExtremeBinDist.Sample(Rnd, newCellsCount - newMutantCount, splitFactor);
                    for (int splitI = 0; splitI < splitCellsCount; splitI++)
                    {
                        var childClone = subClone.CreateChild(GetNewId(), _generation, subClone.DivisionRate, subClone.NumberDrivers);
                        newPops.Add(childClone);
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