using SimChA.Computation;
using SimChA.DataTypes;
using ExtremeBinDist = Extreme.Statistics.Distributions.BinomialDistribution;

namespace SimChA.Simulation;

public class Simulator
{
    private const double MAX_FIT = 10.0;

    public int AliveSC;

    private int newId;

    public int StepNo;

    public Simulator(SimParams simParams, Random rnd)
    {
        StepNo = 0;
        SimParams = simParams;
        Rnd = rnd;

        double initFit = 1;
        double deathRate = 1;
        for (int i = 0; i < SimParams.StartMut; i++)
        {
            double divChange = FitnessFunction.SampleFitness(simParams, rnd);
            initFit = BumpFitness(initFit, SimParams, divChange);
            deathRate = BumpDeath(deathRate, SimParams, divChange);
        }

        var primeval = new SubClone(0, -1, 0, initFit, deathRate, SimParams.StartMut, SimParams.StartPop);
        Clones = new List<SubClone> { primeval };
    }

    public List<SubClone> Clones { get; }
    public SimParams SimParams { get; }
    private Random Rnd { get; }

    private int GetNewId() => ++newId;

    private static double BumpFitness(double original, SimParams simParams, double divChange)
    {
        switch (simParams.FitnessEffect)
        {
            case FitnessEffectType.Death:
                return original;
            case FitnessEffectType.Both:
                divChange /= 2;
                break;
            case FitnessEffectType.Birth:
                break;
        }

        double newFitness = simParams.FitnessAcc switch
        {
            FitnessAccType.Add => original + divChange,
            FitnessAccType.Mul => original * (1 + divChange),
            FitnessAccType.Eth => Math.Clamp(original * (1 + divChange * (1 - original / 10.0)), 0.0, MAX_FIT),
            _ => throw new ArgumentOutOfRangeException()
        };
        return newFitness;
    }

    private static double BumpDeath(double original, SimParams simParams, double divChange)
    {
        switch (simParams.FitnessEffect)
        {
            case FitnessEffectType.Birth:
                return original;
            case FitnessEffectType.Both:
                divChange /= 2;
                break;
            case FitnessEffectType.Death:
                break;
        }

        double newFitness = original / (1 + divChange);
        return newFitness;
    }

    public void Step()
    {
        AliveSC = 0;
        StepNo++;

        List<SubClone> newClones = new();
        long deadCount = CellSampling.NecroCount(Clones);
        long aliveCount = CellSampling.AliveCount(Clones);
        long popSize = deadCount + aliveCount;

        double unconfined = popSize;
        if (SimParams.Confinement > 0)
        {
            double r = Math.Pow(3.0 / 4.0 * (popSize / Math.PI), 1.0 / 3.0);
            double reminder = r - 1.0 / SimParams.Confinement;
            if (reminder > 0)
            {
                double blockedPop = 4.0 / 3.0 * Math.PI * Math.Pow(reminder, 3.0);
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
            int newDead = ExtremeBinDist.Sample(Rnd, (int)subClone.AliveCount, subClone.DeathRate * SimParams.Turnover);
            int newNecrotic = (int)Math.Round(newDead * (1 - divFraction));
            int disappeared = newDead - newNecrotic;

            // Create new cells
            double divRate = Math.Clamp(subClone.BirthRate * divFraction * SimParams.Turnover, 0.0, 1.0);
            int newCellsCount = ExtremeBinDist.Sample(Rnd, (int)subClone.AliveCount, divRate);

            // Mutate some of the cells
            int newMutantCount =
                ExtremeBinDist.Sample(Rnd, newCellsCount, Math.Clamp(SimParams.MutationProb * 2, 0.0, 1.0));

            for (int mutationI = 0; mutationI < newMutantCount; mutationI++)
            {
                double divChange = FitnessFunction.SampleFitness(SimParams, Rnd);
                double newDivision = BumpFitness(subClone.BirthRate, SimParams, divChange);
                double newDeath = BumpDeath(subClone.DeathRate, SimParams, divChange);
                var childClone = subClone.CreateChild(GetNewId(), StepNo, newDivision, newDeath,
                    subClone.NumberDrivers + 1);
                newClones.Add(childClone);
            }

            subClone.NewGen(
                (uint)(subClone.AliveCount + newCellsCount - newMutantCount - newDead),
                (uint)(subClone.NecroCount + newNecrotic),
                (uint)disappeared);
        }

        Clones.AddRange(newClones);
    }
}