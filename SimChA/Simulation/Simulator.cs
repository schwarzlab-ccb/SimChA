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
    

    public void StepTree()
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
    
    public void GetMutations(Clone clone){
        if(Clones.Where(c => c.ParentId == clone.CloneId).FirstOrDefault() != null){
            foreach(var child in Clones.Where(c => c.ParentId == clone.CloneId)){
                child.Karyotype = clone.SetKaryotype();
                var abberation = SelectMutation();
                child.Karyotype.ApplyAbberation(abberation);
                Console.Write("Creating Mutations for Clone " + child.CloneId + "  \r");
                GetMutations(child);
            }
        }
        
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

    public Clone CreateNodes(string newickNode, int parentId){
        string[] cloneString = newickNode.Split(':');
        Clone clone = new Clone(Int32.Parse(cloneString[0].Split('-')[0]), parentId, Int32.Parse(cloneString[1]), Int32.Parse(cloneString[0].Split('-')[1]), new Karyotype(SimParams.IsFemale, Rnd));
        return clone;
    }

    public void BuildCloneFromNewick(string[] newickString){
        Clones.Clear();
        List<int> parentIds = new List<int>();
            parentIds.Add(-1);
            bool rootSet = false;
            for(int i = 0; i < newickString.Count(); i++){
                string test = newickString[i];
                switch(newickString[i]){
                    case "(":
                        if(newickString[i-1] == ""){
                            parentIds = parentIds.Where(p => p != parentIds.Last()).ToList();
                            break;
                        }
                        Clones.Add(CreateNodes(newickString[i-1], parentIds.Last()));
                        parentIds = parentIds.Where(p => p != parentIds.Last()).ToList();
                        break;
                    case ")":
                        if(rootSet){
                            Clones.Add(CreateNodes(newickString[i-1], parentIds.Last()));
                            parentIds.Add(Int32.Parse(newickString[i-1].Split('-')[0]));
                        }
                        else{
                            
                        }
                        break;
                    case ",":
                        if(!rootSet){
                            Clones.Add(CreateNodes(newickString[i-1], parentIds.Last()));
                            parentIds.Add(Int32.Parse(newickString[i-1].Split('-')[0]));
                            rootSet = true;
                        }
                        else
                            Clones.Add(CreateNodes(newickString[i-1], parentIds.Last()));
                        break;
                    default:
                        break;
                }
            }
    }

    public void GetMutationsNewick(Clone newickClone){
        foreach(var clone in Clones.Where(c => c.ParentId==newickClone.CloneId)){
            clone.Karyotype = newickClone.SetKaryotype();
            for(int i = 0; i < clone.MutCount; i++){
                var abberation = SelectMutation();
                clone.Karyotype.ApplyAbberation(abberation);
            }
            GetMutationsNewick(clone);
        }
    }
}