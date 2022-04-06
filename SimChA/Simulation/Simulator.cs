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
        double initFit = SimParams.BirthRate;
        for (int i = 0; i < SimParams.InitMut; i++)
        {
            initFit = BumpFitness(initFit, SimParams, Rnd);
        }

        var primeval = new SubClone(0, -1, 0, initFit);
        // var primeval = new SubClone(0, -1, 0, BumpFitness(SimParams.BirthRate, SimParams, Rnd));
        Clones = new List<SubClone> {primeval};
    }

    private static double BumpFitness(double original, SimParams simParams, Random rnd)
    {
        double divChange = FitnessFunction.SampleFitness(simParams, rnd);
        double newFitness = simParams.MultiplicativeFitness
            ? original * (1 + divChange)
            : original + divChange * simParams.BirthRate;
        return newFitness;
    }

    public void Step()
    {
        _generation++;
        AliveSC = 0;

        List<SubClone> newClones = new();
        long popSize = CellSampling.PopulationSize(Clones);
        long aliveCount = CellSampling.AliveCount(Clones);
        double divisionFraction = 1;
        if (SimParams.Confinement > 0 && aliveCount > 1 / SimParams.Confinement)
        {
            double divisible = Math.Round(Math.Pow(popSize, TWO_THIRD)) / SimParams.Confinement;
            if (aliveCount > divisible && aliveCount > 0)
            {
                divisionFraction = Math.Clamp(divisible / aliveCount, 0.0, 1.0);
            }
        }

        foreach (var subClone in Clones.Where(sc => sc.AliveCount > 0))
        {
            AliveSC++;

            // Kill cells
            int newDead;
            if (SimParams.StochasticCellLife)
            {
                // newDead = Binomial.Sample(Rnd, SimParams.DivisionRate, (int)subClone.AliveCount);
                newDead = ExtremeBinDist.Sample(Rnd, (int) subClone.AliveCount, SimParams.BirthRate);
            }
            else
            {
                double deadFraction = SimParams.BirthRate * subClone.AliveCount + subClone.ToDie;
                newDead = (int) deadFraction;
                subClone.ToDie = deadFraction - newDead;
            }

            // if (subClone.TotalCount <= SimParams.InitPop)
            // {
            //     newDead = Math.Min(newDead, (int)subClone.AliveCount - 1);
            // }

            // Create new cells
            int newCellsCount;
            double divRate = Math.Clamp(subClone.DivisionRate * divisionFraction, 0.0, 1.0);
            if (SimParams.StochasticCellLife)
            {
                newCellsCount =
                    ExtremeBinDist.Sample(Rnd, (int) subClone.AliveCount, divRate);
                // newCellsCount = Binomial.Sample(Rnd, divRate, (int)subClone.AliveCount);
            }
            else
            {
                double divideFraction = (divRate > 0 ? divRate * subClone.AliveCount : 0) + subClone.ToDivide;
                newCellsCount = (int) divideFraction;
                subClone.ToDivide = divideFraction - newCellsCount;
            }

            // Mutate some of the cells
            int newMutantCount = SimParams.MutationRate > 0
                ? ExtremeBinDist.Sample(Rnd, newCellsCount * 2, SimParams.MutationRate)
                : 0;

            for (int mutationI = 0; mutationI < newMutantCount; mutationI++)
            {
                double newDivision = BumpFitness(subClone.DivisionRate, SimParams, Rnd);
                var childClone = subClone.CreateChild(GetNewId(), _generation, newDivision, subClone.NumberDrivers + 1);
                newClones.Add(childClone);
            }


            subClone.NewGen(
                (uint) (subClone.AliveCount + newCellsCount - newMutantCount - newDead),
                (uint) (subClone.DeadCount + newDead));
        }

        Clones.AddRange(newClones);
    }
}