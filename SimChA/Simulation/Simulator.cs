using Extreme.Statistics.Distributions;
using SimChA.DataTypes;
using ExtremeBinDist = Extreme.Statistics.Distributions.BinomialDistribution;

namespace SimChA.Simulation;

public static class Simulator
{
    public static void GetMutationsNewick(Clone newickClone, List<Clone> clones, AberrationsInfo aberrationsInfo, Random rnd)
    {
        // TODO this is still square complexity, fix! (should use a tree structure)
        foreach (var clone in clones.Where(c => c.ParentId == newickClone.CloneId))
        {
            clone.Karyotype = newickClone.SetKaryotype();
            for (int i = 0; i < clone.MutCount; i++)
            {
                var aberration = aberrationsInfo.PickRandomMutation(rnd);
                clone.Karyotype.ApplyAberration(rnd, aberration, aberrationsInfo.Map[aberration]);
            }
            GetMutationsNewick(clone, clones, aberrationsInfo, rnd);
        }
    }
}