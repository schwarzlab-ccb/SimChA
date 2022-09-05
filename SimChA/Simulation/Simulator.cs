using Extreme.Statistics.Distributions;
using SimChA.DataTypes;
using SimChA.IO;
using ExtremeBinDist = Extreme.Statistics.Distributions.BinomialDistribution;

namespace SimChA.Simulation;

public class Simulator
{
    public Simulator(SimParams simParams, Random rnd)
    {
        Rnd = rnd;
        LastId = -1;
        Aberrations = new Aberrations(simParams);
        IsFemale = simParams.IsFemale;
        Clones = new List<Clone>();
    }

    private bool IsFemale;
    public List<Clone> Clones { get; }
    public Aberrations Aberrations { get; }
    private Random Rnd { get; }
    private int LastId { get; set; }
    private int GetNewId() => ++LastId;

    private AberrationEnum SelectMutation()
    {
        double ratesSum = Aberrations.RatesSum;
        double sample = ContinuousUniformDistribution.Sample(Rnd, 0, ratesSum);
        foreach (var rate in Aberrations.Map)
        {
            if (sample <= rate.Value.Likelihood)
            {
                return rate.Key;
            }
            sample -= rate.Value.Likelihood;
        }

        // In case float-point calculations would cause jumping out of the loop
        return Aberrations.Map.Last().Key;
    }

    private Clone CreateNodes(string newickNode, int parentId)
    {
        string[] cloneString = newickNode.Split(':');
        // TODO: split below in individual assignments
        var clone = new Clone(int.Parse(cloneString[0].Split('-')[0]), parentId, int.Parse(cloneString[1]),
            int.Parse(cloneString[0].Split('-')[1]), new Karyotype(IsFemale, Rnd));
        return clone;
    }

    public void BuildCloneFromNewick(string[] newickString)
    {
        var parentIds = new List<int> { -1 };
        bool rootSet = false;
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

                    Clones.Add(CreateNodes(newickString[i - 1], parentIds.Last()));
                    parentIds = parentIds.Where(p => p != parentIds.Last()).ToList();
                    break;
                case ")":
                    if (rootSet)
                    {
                        Clones.Add(CreateNodes(newickString[i - 1], parentIds.Last()));
                        parentIds.Add(int.Parse(newickString[i - 1].Split('-')[0]));
                    }

                    break;
                case ",":
                    if (!rootSet)
                    {
                        Clones.Add(CreateNodes(newickString[i - 1], parentIds.Last()));
                        parentIds.Add(int.Parse(newickString[i - 1].Split('-')[0]));
                        rootSet = true;
                    }
                    else
                    {
                        Clones.Add(CreateNodes(newickString[i - 1], parentIds.Last()));
                    }

                    break;
            }
        }
    }

    public void GetMutationsNewick(Clone newickClone)
    {
        foreach (var clone in Clones.Where(c => c.ParentId == newickClone.CloneId))
        {
            clone.Karyotype = newickClone.SetKaryotype();
            for (int i = 0; i < clone.MutCount; i++)
            {
                var aberration = SelectMutation();
                clone.Karyotype.ApplyAbberation(aberration);
            }

            GetMutationsNewick(clone);
        }
    }
}