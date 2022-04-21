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

    public Simulator(SimParams simParams, Random rnd)
    {
        StepNo = 0;
        SimParams = simParams;
        Rnd = rnd;

        double initFit = 1;
        for (int i = 0; i < SimParams.InitMut; i++)
        {
            initFit = BumpFitness(initFit, SimParams, Rnd);
        }

        var primeval = new SubClone(0, -1, 0, initFit, SimParams.InitMut, SimParams.InitPop);
        Clones = new List<SubClone> { primeval };
    }

    private static double BumpFitness(double original, SimParams simParams, Random rnd)
    {
        double divChange = FitnessFunction.SampleFitness(simParams, rnd);
        double newFitness = simParams.FitnessAcc switch
        {
            FitnessAccType.Add => original + divChange,
            FitnessAccType.Mul => original * (1 + divChange),
            FitnessAccType.Eth => original * (1 + divChange * ( 1 - original / 10.0)),
            _ => throw new ArgumentOutOfRangeException()
        };
        return newFitness;
    }

    public void Step()
    {
        AliveSC = 0;
        StepNo++;

        List<SubClone> newClones = new();
        long deadCount = CellSampling.DeadCount(Clones);
        long aliveCount = CellSampling.AliveCount(Clones);
        long popSize = deadCount + aliveCount;
        
        double unconfined = popSize;
        if (SimParams.Confinement > 0)
        {
            double r = Math.Pow((3.0 * popSize) / (4.0 * Math.PI), 1.0 / 3.0);
            double reminder = r - (1 / SimParams.Confinement);
            if (reminder > 0)
            {
                double blockedPop = 4.0 / 3.0 * Math.PI * Math.Pow(reminder, 3);
                unconfined = popSize - blockedPop;
            }
        }
        
        double divFraction = aliveCount > unconfined && aliveCount > 0
            ? Math.Clamp(unconfined / aliveCount, 0.0, 1.0)
            : 1.0;

        foreach (var subClone in Clones.Where(sc => sc.AliveCount > 0))
        {
            AliveSC++;

            // Kill cells
            int newDead = ExtremeBinDist.Sample(Rnd, (int)subClone.AliveCount, SimParams.Turnover);
            int newNecrotic = (int) Math.Round(newDead * (1 - divFraction));
            int disappeared = newDead - newNecrotic;

            // Create new cells
            double divRate = Math.Clamp(subClone.DivisionRate * divFraction * SimParams.Turnover, 0.0, 1.0);
            int newCellsCount = ExtremeBinDist.Sample(Rnd, (int)subClone.AliveCount, divRate);

            // Mutate some of the cells
            int newMutantCount = ExtremeBinDist.Sample(Rnd, newCellsCount, Math.Clamp(SimParams.MutationProb * 2, 0.0, 1.0));

            for (int mutationI = 0; mutationI < newMutantCount; mutationI++)
            {
                double newDivision = BumpFitness(subClone.DivisionRate, SimParams, Rnd);
                var childClone = subClone.CreateChild(GetNewId(), StepNo, newDivision, subClone.NumberDrivers + 1);
                newClones.Add(childClone);
            }

            subClone.NewGen(
                (uint)(subClone.AliveCount + newCellsCount - newMutantCount - newDead),
                (uint)(subClone.DeadCount + newNecrotic),
                (uint) disappeared);
        }

        Clones.AddRange(newClones);
    }
}