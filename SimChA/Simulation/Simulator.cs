using MathNet.Numerics.Distributions;
using SimChA.Computation;
using SimChA.DataTypes;
using ExtremeBinDist = Extreme.Statistics.Distributions.BinomialDistribution;

namespace SimChA.Simulation;

public class Simulator
{
    public List<SubClone> Clones { get; }
    public SimParams SimParams { get; }

    public int AliveSC;

    private int newId;
    private int GetNewId() => ++newId;

    public int StepNo;
    private Random Rnd { get; }

    private const int INIT_POP = 1;
    const double THIRD = 1.0 / 3.0;
    const double TWO_THIRD = 2.0 / 3.0;

    public Simulator(SimParams simParams, Random rnd)
    {
        StepNo = 0;
        SimParams = simParams;
        Rnd = rnd;

        double initFit = SimParams.Turnover;
        for (int i = 0; i < SimParams.InitMut; i++)
        {
            initFit = BumpFitness(initFit, SimParams, Rnd);
        }
        var primeval = new SubClone(0, -1, 0, initFit);
        Clones = new List<SubClone> {primeval};
    }

    private static double BumpFitness(double original, SimParams simParams, Random rnd)
    {
        double divChange = FitnessFunction.SampleFitness(simParams, rnd);
        double newFitness = simParams.MultiplicativeFitness
            ? original * (1 + divChange)
            : original + divChange * simParams.Turnover;
        return newFitness;
    }

    public void Step()
    {
        AliveSC = 0;
        StepNo++;

        List<SubClone> newClones = new();
        long popSize = CellSampling.PopulationSize(Clones);
        long aliveCount = CellSampling.AliveCount(Clones);
        double divisionFraction = 1;
        if (SimParams.Confinement > 0)
        {
            double divisible = Math.Round(Math.Pow(popSize, TWO_THIRD)) / SimParams.Confinement;
            if (aliveCount > divisible && aliveCount > 0)
            {
                divisionFraction = Math.Clamp(divisible / aliveCount, 0.0, 1.0);
            }
        }

        foreach (var subClone in Clones.Where(sc => sc.SampleCount > 0))
        {
            AliveSC++;

            // Kill cells
            double deathRate = Math.Clamp(SimParams.Turnover * (1 - divisionFraction), 0.0, 1.0);
            int newDead = ExtremeBinDist.Sample(Rnd, (int) subClone.AliveCount, deathRate);
            
            // Create new cells
            double divRate = Math.Clamp(subClone.DivisionRate * divisionFraction, 0.0, 1.0);
            int newCellsCount = ExtremeBinDist.Sample(Rnd, (int) subClone.AliveCount, divRate);

            // Decay cells 
            int decayedCount = ExtremeBinDist.Sample(Rnd, (int) subClone.DeadCount, SimParams.Turnover);

            // Mutate some of the cells
            int newMutantCount = ExtremeBinDist.Sample(Rnd, newCellsCount * 2, SimParams.MutationProb);

            for (int mutationI = 0; mutationI < newMutantCount; mutationI++)
            {
                double newDivision = BumpFitness(subClone.DivisionRate, SimParams, Rnd);
                var childClone = subClone.CreateChild(GetNewId(), StepNo, newDivision, subClone.NumberDrivers + 1);
                newClones.Add(childClone);
            }

            subClone.NewGen(
                (uint) (subClone.AliveCount + newCellsCount - newMutantCount - newDead),
                (uint) (subClone.DeadCount + newDead - decayedCount),
                (uint) (subClone.DecayedCount + decayedCount));
        }

        Clones.AddRange(newClones);
    }
}