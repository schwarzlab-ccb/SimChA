using MathNet.Numerics.Distributions;
using SimChA.DataTypes;

namespace SimChA.Simulation;

public class Simulator
{
    public List<SubClone> Clones { get; }
    public SimParams SimParams { get; }

    public Simulator(SimParams simParams)
    {
        var firstClone = new SubClone(0, -1, simParams.DivisionRate, simParams.MutationRate,
            new Karyotype(simParams.IsFemale));
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
            int newMutantCount =
                Binomial.Sample(originalClone.MutationRate, newCellsCount); // The existing cells will not mutate
            for (int mutationI = 0; mutationI < newMutantCount; mutationI++)
            {
                var newSubClone = new SubClone(originalClone, Clones.Count);
                var abberation = SelectMutation();
                ApplyMutation(newSubClone.Karyotype, abberation);
                Clones.Add(newSubClone);
            }

            originalClone.AliveCount += newCellsCount - newMutantCount;
        }
    }

    private AbberationEnum SelectMutation()
    {
        return AbberationEnum.TailDeletion;
    }

    private void ApplyMutation(Karyotype karyotype, AbberationEnum abberation)
    {
        switch (abberation)
        {
            case AbberationEnum.TailDeletion:
                karyotype.ApplyTailDeletion();
                break;
            case AbberationEnum.Missegregation:
            case AbberationEnum.Duplication:
            case AbberationEnum.Chromothripsis:
            case AbberationEnum.Translocation:
            case AbberationEnum.InternalDeletion:
            case AbberationEnum.Inversion:
            case AbberationEnum.BreakageFusionBridge:
            default:
                throw new ArgumentOutOfRangeException();
        }
    }
}