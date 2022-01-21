using MathNet.Numerics.Distributions;
using SimChA.DataTypes;

namespace SimChA.Simulation;

public class Simulator
{
    public List<SubClone> Clones { get; }
    public SimParams SimParams { get; }

    private int _generation;

    public Simulator(SimParams simParams)
    {
        _generation = 0;
        var refKaryotype = new Karyotype(simParams.IsFemale);
        var firstClone = new SubClone(0, -1, _generation, simParams.DivisionRate, simParams.MutationRate, refKaryotype, simParams.InitialPop);
        Clones = new List<SubClone> { firstClone };
        SimParams = simParams;
    }

    public void Step()
    {
        _generation++;
        Kill();
        DivideAndMutate();
    }

    private void Kill()
    {
        foreach (var subClone in Clones.Where(sc => sc.AliveCount > 0))
        {
            int currentPop = subClone.Generations[_generation - 1];
            int deadCount = Binomial.Sample(SimParams.DeathRate, currentPop);
            subClone.Generations[_generation] = currentPop - deadCount;
        }
    }

    private bool IsViable(Karyotype kar) => kar.ChromCount > 23;

    private void DivideAndMutate()
    {
        List<SubClone> newClones = new();
        foreach (var subClone in Clones.Where(sc => sc.AliveCount > 0))
        {
            int newCellsCount = Binomial.Sample(subClone.DivisionRate, subClone.AliveCount);
            // The existing cells will not mutate
            int newMutantCount = Binomial.Sample(subClone.MutationRate, newCellsCount); 
            for (int mutationI = 0; mutationI < newMutantCount; mutationI++)
            {
                var childClone = subClone.CreateChild(Clones.Count + newClones.Count, _generation);
                var abberation = SelectMutation();
                childClone.Karyotype.ApplyAbberation(abberation);
                if (!IsViable(childClone.Karyotype)) 
                    continue;
                
                childClone.DivisionRate = Math.Clamp(childClone.DivisionRate * SimParams.FitnessInc, 0, 1);
                newClones.Add(childClone);
            }
            subClone.Generations[_generation] += newCellsCount - newMutantCount;
        }
        Clones.AddRange(newClones);
    }

    private AbberationEnum SelectMutation()
    {
        double ratesSum = SimParams.RatesSum;
        double sample = ContinuousUniform.Sample(0, ratesSum);
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