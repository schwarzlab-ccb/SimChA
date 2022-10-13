using Extreme.Statistics.Distributions;
using SimChA.DataTypes;
using ExtremeBinDist = Extreme.Statistics.Distributions.BinomialDistribution;

namespace SimChA.Simulation;

public static class Simulator
{
    public static void AssignMutationsRecursive(Clone currentClone, List<Clone> clones, AberrationsInfo aberrationsInfo, Random rnd)
    {
        // TODO this is still square complexity, fix! (should use a tree structure)
        foreach (var clone in clones.Where(c => c.ParentId == currentClone.CloneId))
        {
            clone.Karyotype = currentClone.SetKaryotype();
            for (int i = 0; i < clone.MutCount; i++)
            {
                var aberration = aberrationsInfo.PickRandomMutation(rnd);
                clone.Karyotype.ApplyAberration(rnd, aberration, aberrationsInfo.Map[aberration]);
            }
            AssignMutationsRecursive(clone, clones, aberrationsInfo, rnd);
        }
    }

    public static List<Clone> GetClonePair(int distance, bool isFemale)
    {
        var parent = new Clone(0, -1, "1", 0, 0, new Karyotype(isFemale));
        var child = new Clone(1, 0, "2", distance, 1, new Karyotype(isFemale));
        return new List<Clone> {parent, child};
    }
}