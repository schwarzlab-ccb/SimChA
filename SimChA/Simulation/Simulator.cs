using SimChA.Computation;
using SimChA.DataTypes;
using static Extreme.Statistics.Distributions.BinomialDistribution;
using ExtremeBinDist = Extreme.Statistics.Distributions.BinomialDistribution;

namespace SimChA.Simulation;

public class Simulator
{
    public List<Clone> Clones { get; }
    public SimParams SimParams { get; }
    private Random Rnd { get; }
    public int AliveClones => Clones.Count(c => c.IsAlive);
    public int StepNo { get; private set; }
    private int LastId { get; set; }
    private int GetNewId() => ++LastId;
    public long CellCount => Clones.Sum(clone => clone.CellCount);

    public Simulator(SimParams simParams, Random rnd)
    {
        Rnd = rnd;
        StepNo = 0;
        LastId = -1;
        SimParams = simParams;
        var initialKaryotype = new Karyotype(simParams.IsFemale, rnd);
        var primeval = new Clone(GetNewId(), -1, SimParams.StartMut, SimParams.StartPop, initialKaryotype);
        Clones = new List<Clone> { primeval };
    }
    
    public void Step()
    {
        StepNo++;

        List<Clone> newClones = new();

        foreach (var clone in Clones.Where(c => c.IsAlive))
        {
            int newDead = Sample(Rnd, clone.CellCount, SimParams.DeathRate * SimParams.Turnover);
            int divisions = Sample(Rnd, clone.CellCount, SimParams.Turnover);
            int mutations = Sample(Rnd, 2 * divisions, SimParams.MutationProb);
            clone.CellCount = Math.Max(0, clone.CellCount + divisions - newDead - mutations);
            for (int i = 0; i < mutations; i++)
            {
                var newClone = clone.CreateChild(GetNewId());
                var abberation = SelectMutation();
                newClone.Karyotype.ApplyAbberation(abberation);
                newClones.Add(newClone);
            }
        }

        Clones.AddRange(newClones);
    }
    
    private AberrationEnum SelectMutation()
    {
        double ratesSum = SimParams.SumRates();
        double sample = Extreme.Statistics.Distributions.ContinuousUniformDistribution.Sample(Rnd, 0, ratesSum);
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
}