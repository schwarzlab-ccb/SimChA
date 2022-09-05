using Extreme.Statistics.Distributions;
using SimChA.DataTypes;
using ExtremeBinDist = Extreme.Statistics.Distributions.BinomialDistribution;

namespace SimChA.Simulation;

public static class Simulator
{
    private static AberrationEnum SelectMutation(AberrationsInfo aberrationsInfo, Random rnd)
    {
        double ratesSum = aberrationsInfo.RatesSum;
        double sample = ContinuousUniformDistribution.Sample(rnd, 0, ratesSum);
        foreach (var rate in aberrationsInfo.Map)
        {
            if (sample <= rate.Value.Likelihood)
            {
                return rate.Key;
            }
            sample -= rate.Value.Likelihood;
        }

        // In the case that float-point calculations would cause jumping out of the loop, use the last one
        return aberrationsInfo.Map.Last().Key;
    }

    private static Clone CreateNodes(string newickNode, int parentId, bool isFemale, Random rnd)
    {
        string[] cloneString = newickNode.Split(':');
        // TODO: split below in individual assignments
        var clone = new Clone(int.Parse(cloneString[0].Split('-')[0]), parentId, int.Parse(cloneString[1]),
            int.Parse(cloneString[0].Split('-')[1]), new Karyotype(isFemale, rnd));
        return clone;
    }

    public static List<Clone> BuildCloneFromNewick(string[] newickString, bool isFemale, Random rnd)
    {
        List<Clone> clones = new();
        var parentIds = new List<int> { -1 };
        bool rootSet = false;
        // TODO: Multiple code repetitions below, fix
        for (int i = 0; i < newickString.Length; i++)
        {
            switch (newickString[i])
            {
                case "(":
                    if (newickString[i - 1] == "")
                    {
                        parentIds = parentIds.Where(p => p != parentIds.Last()).ToList();
                        break;
                    }

                    clones.Add(CreateNodes(newickString[i - 1], parentIds.Last(), isFemale, rnd));
                    parentIds = parentIds.Where(p => p != parentIds.Last()).ToList();
                    break;
                case ")":
                    if (rootSet)
                    {
                        clones.Add(CreateNodes(newickString[i - 1], parentIds.Last(), isFemale, rnd));
                        parentIds.Add(int.Parse(newickString[i - 1].Split('-')[0]));
                    }

                    break;
                case ",":
                    if (!rootSet)
                    {
                        clones.Add(CreateNodes(newickString[i - 1], parentIds.Last(), isFemale, rnd));
                        parentIds.Add(int.Parse(newickString[i - 1].Split('-')[0]));
                        rootSet = true;
                    }
                    else
                    {
                        clones.Add(CreateNodes(newickString[i - 1], parentIds.Last(), isFemale, rnd));
                    }

                    break;
            }
        }
        return clones;
    }

    public static void GetMutationsNewick(Clone newickClone, List<Clone> clones, AberrationsInfo aberrationsInfo, Random rnd)
    {
        foreach (var clone in clones.Where(c => c.ParentId == newickClone.CloneId))
        {
            clone.Karyotype = newickClone.SetKaryotype();
            for (int i = 0; i < clone.MutCount; i++)
            {
                var aberration = SelectMutation(aberrationsInfo, rnd);
                clone.Karyotype.ApplyAberration(aberration, aberrationsInfo.Map[aberration]);
            }
            GetMutationsNewick(clone, clones, aberrationsInfo, rnd);
        }
    }
}