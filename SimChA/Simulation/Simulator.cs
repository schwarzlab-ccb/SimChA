using Extreme.Statistics.Distributions;
using SimChA.DataTypes;
using ExtremeBinDist = Extreme.Statistics.Distributions.BinomialDistribution;

namespace SimChA.Simulation;

public static class Simulator
{
    public static void AssignMutationsRecursive(Clone currentClone, List<Clone> clones, AberrationsInfo aberrationsInfo, Random rnd, List<Gen> essential, List<Gen> tsgOg)
    {
        // TODO this is still square complexity, fix! (should use a tree structure)
        foreach (int cloneId in currentClone.ChildrenIDs)
        {
            clones[cloneId].Karyotype = currentClone.SetKaryotype();
            for (int i = 0; i < clones[cloneId].MutCount; i++)
            {
                var aberration = aberrationsInfo.PickRandomMutation(rnd);
                clones[cloneId].deltaFitness = clones[cloneId].Karyotype.ApplyAberration(rnd, aberration, aberrationsInfo.Map[aberration], clones[cloneId].Name);
            }
            AssignFitness(clones[cloneId], essential, tsgOg);
            AssignMutationsRecursive(clones[cloneId], clones, aberrationsInfo, rnd, essential, tsgOg);
        }
    }

    public static List<Clone> GetClonePair(int distance, bool isFemale)
    {
        var parent = new Clone(0, -1, "1", 0, new Karyotype(isFemale));
        var child = new Clone(1, 0, "2", distance, new Karyotype(isFemale));
        parent.ChildrenIDs.Add(child.CloneId);
        return new List<Clone> {parent, child};
    }

    public static void AssignFitness(Clone clone, List<Gen> essentialList, List<Gen> tsgOgList)
    {
        clone.deltaFitness = 0;
        float stress = 0;
        float essentialityFitness = 0;
        float tsgOgFitness = 0;
        if(clone.Karyotype.ChromCount > 46)
        {
            stress = (float)(-0.0001 * Math.Pow(clone.Karyotype.ChromCount-46, 2));
        }

        List<Region> regionsLost = new List<Region>();
        List<Gen> essentialStillThere = new List<Gen>();
        List<Gen> tsgOgStillThere = new List<Gen>();
        foreach(Chromosome chr in clone.Karyotype.Chromosomes)
        {
            for(int i = 0; i < chr._regions.Count(); i++)
            {
                essentialStillThere.AddRange(essentialList.FindAll(x => x.start > chr._regions[i].Start && x.stop < chr._regions[i].End && x.chr == chr._regions[i].ChromId.ChromNum));
                tsgOgStillThere.AddRange(tsgOgList.FindAll(x => x.start > chr._regions[i].Start && x.stop < chr._regions[i].End && x.chr == chr._regions[i].ChromId.ChromNum));
            }
        }     
        List<Gen> essentialsLost = essentialList.Except(essentialStillThere).ToList();
        List<Gen> tsgOgLost = tsgOgList.Except(tsgOgStillThere).ToList();
        List<Gen> tsgOgAdded = new List<Gen>();
        foreach(var tsgOg in tsgOgStillThere.GroupBy(x => x.name))
        {
            for(int i = 2; i < tsgOg.Count(); i++){
                tsgOgAdded.Add(tsgOg.ElementAt(0));
            }
        }
        clone.deltaFitness = stress;
    }
}