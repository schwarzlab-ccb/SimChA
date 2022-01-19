using MathNet.Numerics.Distributions;
using SimChA.DataTypes;

namespace SimChA.Simulation;

public class Simulator
{
    public List<SubClone> Clones { get; }
    public SimParams SimParams { get; }

    public Simulator(SimParams simParams)
    {
        var initKaryotype = new Karyotype(simParams.IsFemale);
        var firstClone = new SubClone(0, -1, simParams.DivisionRate, simParams.MutationRate, initKaryotype);
        Clones = new List<SubClone> { firstClone };
        SimParams = simParams;
    }

    public void Step()
    {
        DivideAndMutate();
        // Kill();
    }

    private void Kill()
    {
        throw new NotImplementedException();
    }

    private void DivideAndMutate()
    {
        int cloneCount = Clones.Count;
        for (int i = 0; i < cloneCount; i++)
        {
            var originalClone = Clones[i];
            int newCellsCount = Binomial.Sample(originalClone.DivisionRate, originalClone.AliveCount);
            // The existing cells will not mutate
            int newMutantCount = Binomial.Sample(originalClone.MutationRate, newCellsCount); 
            for (int mutationI = 0; mutationI < newMutantCount; mutationI++)
            {
                var newSubClone = originalClone.CreateChild(Clones.Count);
                var abberation = SelectMutation();
                newSubClone.Karyotype.ApplyAbberation(abberation);
                if (newSubClone.Karyotype.ChromCount > 23)
                {
                    Clones.Add(newSubClone);
                }
            }

            originalClone.AliveCount += newCellsCount - newMutantCount;
        }
    }

    private AbberationEnum SelectMutation()
    {
        double ratesSum = SimParams.RatesSum;
        double sample = ContinuousUniform.Sample(0, ratesSum);
        foreach (var rate in SimParams.AbberationRates)
        {
            if (sample <= rate.Value) 
            {
                return rate.Key;
            }
            sample -= rate.Value;
        }
        // In case float-point calculations would cause jumping out of the loop
        return SimParams.AbberationRates.Last().Key;
    }
}