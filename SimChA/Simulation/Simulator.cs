using Extreme.Statistics.Distributions;
using SimChA.DataTypes;
using ExtremeBinDist = Extreme.Statistics.Distributions.BinomialDistribution;

namespace SimChA.Simulation;

public class Simulator
{
    public Simulator(SimParams simParams, Random rnd)
    {
        Rnd = rnd;
        LastId = -1;
        SimParams = simParams;
        var initialKaryotype = new Karyotype(simParams.IsFemale, rnd);
        var primeval = new Clone(GetNewId(), -1, 0, 0, initialKaryotype);
        Clones = new List<Clone> { primeval };
    }

    public List<Clone> Clones { get; }
    public SimParams SimParams { get; }
    private Random Rnd { get; }
    private int LastId { get; set; }
    private int GetNewId() => ++LastId;


    private AberrationEnum SelectMutation()
    {
        double ratesSum = SimParams.SumRates();
        double sample = ContinuousUniformDistribution.Sample(Rnd, 0, ratesSum);
        foreach (var rate in SimParams.AberrationRates)
        {
            if (sample <= rate.Value)
            {
                return rate.Key;
            }

            sample -= rate.Value;
        }

        // In case float-point calculations would cause jumping out of the loop
        return SimParams.AberrationRates.Last().Key;
    }

    private Clone CreateNodes(string newickNode, int parentId)
    {
        string[] cloneString = newickNode.Split(':');
        var clone = new Clone(int.Parse(cloneString[0].Split('-')[0]), parentId, int.Parse(cloneString[1]),
            int.Parse(cloneString[0].Split('-')[1]), new Karyotype(SimParams.IsFemale, Rnd));
        return clone;
    }

    public void BuildCloneFromNewick(string[] newickString)
    {
        Clones.Clear();
        var parentIds = new List<int>();
        parentIds.Add(-1);
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