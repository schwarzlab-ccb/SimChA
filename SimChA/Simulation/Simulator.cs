using Extreme.Statistics.Distributions;
using SimChA.DataTypes;
using ExtremeBinDist = Extreme.Statistics.Distributions.BinomialDistribution;
using SimChA.IO;

namespace SimChA.Simulation;

public static class Simulator
{
    public static void AssignMutationsRecursive(Clone currentClone, List<Clone> clones, AberrationsInfo aberrationsInfo, Random rnd)
    {
        foreach (int cloneId in currentClone.ChildrenIDs)
        {
            clones[cloneId].Karyotype = currentClone.SetKaryotype();
            for (int i = 0; i < clones[cloneId].MutCount; i++)
            {
                var aberration = aberrationsInfo.PickRandomMutation(rnd);
                string region = clones[cloneId].Karyotype.ApplyAberration(rnd, aberration, aberrationsInfo.Map[aberration], clones[cloneId].Name);
                AssignFitness(clones[cloneId]);
                new Abberation(clones[cloneId].Name, aberration.ToString(), region, (float) clones[cloneId].deltaFitness);
            }
            AssignMutationsRecursive(clones[cloneId], clones, aberrationsInfo, rnd);
        }
    }

    public static List<Clone> GetClonePair(int distance, bool isFemale)
    {
        var parent = new Clone(0, -1, "1", 0, new Karyotype(isFemale));
        var child = new Clone(1, 0, "2", distance, new Karyotype(isFemale));
        parent.ChildrenIDs.Add(child.CloneId);
        return new List<Clone> {parent, child};
    }

    public static void AssignFitness(Clone clone)
    {
        clone.deltaFitness = 0;
        float stress = 0;
        float essentialityFitness = 0;
        float tsgOgFitness = 0;
        if(clone.Karyotype.ChromCount > 46)
        {
            stress = (float)(-0.0001 * Math.Pow(clone.Karyotype.ChromCount-46, 2));
        }
        List<Gen> essentialStillThere = new List<Gen>();
        List<Gen> tsgOgStillThere = new List<Gen>();
        foreach(Chromosome chr in clone.Karyotype.Chromosomes)
        {
            for(int i = 0; i < chr._regions.Count(); i++)
            {

                essentialStillThere.AddRange(GenList.EssentialList[chr._regions[i].ChromId.ChromNum].FindAll(x => x.start > chr._regions[i].Start && x.stop < chr._regions[i].End));
                tsgOgStillThere.AddRange(GenList.TsgOgList[chr._regions[i].ChromId.ChromNum].FindAll(x => x.start > chr._regions[i].Start && x.stop < chr._regions[i].End));
            }
        }     

        List<Gen> essentialsLost = FindMissingGenes(essentialStillThere, GenList.EssentialList);
        List<Gen> tsgOgLost = FindMissingGenes(tsgOgStillThere, GenList.TsgOgList);
        tsgOgLost.AddRange(tsgOgStillThere.GroupBy(x => x).Where(g => g.Count() < 2).Select(g => g.Key));
        List<Gen> tsgOgAdded = new List<Gen>();
        tsgOgAdded.AddRange(tsgOgStillThere.GroupBy(x => x).Where(g => g.Count() > 2).Select(global => global.Key));
        foreach(var essential in essentialsLost)
        {
            essentialityFitness += essential.deltaFitness;
        }
        foreach(var tsgOg in tsgOgAdded)
        {
            tsgOgFitness -= tsgOg.deltaFitness;
        }
        foreach(var tsgOg in tsgOgLost)
        {
            tsgOgFitness += tsgOg.deltaFitness;
        }
        
        clone.deltaFitness = stress + tsgOgFitness + essentialityFitness;
    }

    public static bool AddIfNotThere<ChromNum, Gen>(this Dictionary<ChromNum, List<Gen>> dict, ChromNum key, Gen value)
    {
        if(dict.ContainsKey(key))
        {
            dict[key].Add(value);
        }
        else
        {
            dict.Add(key, new List<Gen> {value});
        }
        return true;
    }

    private static List<Gen> FindMissingGenes(List<Gen> genes, Dictionary<ChromNum, List<Gen>> dict)
    {
        List<Gen> missingGenes = new List<Gen>();
        foreach(var list in dict){
            missingGenes.AddRange(list.Value.Except(genes));
        }
        return missingGenes;
    }
}