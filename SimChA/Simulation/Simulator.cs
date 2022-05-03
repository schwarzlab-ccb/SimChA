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
            initFit = BumpFitness(initFit, SimParams, divChange, true);
            deathRate = BumpFitness(deathRate, SimParams, divChange, false);
        }

        var primeval = new SubClone(0, -1, 0, initFit, deathRate, SimParams.StartMut, SimParams.StartPop);
        Clones = new List<SubClone> { primeval };
    }

    public List<SubClone> Clones { get; }
    public SimParams SimParams { get; }
    private Random Rnd { get; }

    private int GetNewId() => ++newId;

    private static double AccFitness(double original, double change, FitnessAccType type) 
        => type switch
        {
            FitnessAccType.Add => original + change,
            FitnessAccType.Mul => original * (1 + change),
            FitnessAccType.Eth => Math.Clamp(original * (1 + change * (1 - original / MAX_FIT)), 0.0, MAX_FIT),
            _ => original
        };

    private static double BumpFitness(double original, SimParams simParams, double divChange, bool isBirth) 
        => simParams.FitnessEffect switch
        {
            FitnessEffectType.Birth when isBirth => AccFitness(original, divChange, simParams.FitnessAcc),
            FitnessEffectType.Death when !isBirth => AccFitness(original, divChange, simParams.FitnessAcc),
            FitnessEffectType.Both => AccFitness(original, divChange / 2, simParams.FitnessAcc),
            FitnessEffectType.Death when isBirth => original,
            FitnessEffectType.Birth when !isBirth => original,
            _ => original
        };

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
            int newDead = ExtremeBinDist.Sample(Rnd, (int)subClone.AliveCount, (1 / subClone.DeathRate) * SimParams.Turnover);
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
                double newDivision = BumpFitness(subClone.BirthRate, SimParams, divChange, true);
                double newDeath = BumpFitness(subClone.DeathRate, SimParams, divChange, false);
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