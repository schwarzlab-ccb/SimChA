using Extreme.Statistics.Distributions;
using SimChA.DataTypes;
using ExtremeBinDist = Extreme.Statistics.Distributions.BinomialDistribution;
using SimChA.IO;

namespace SimChA.Simulation;

public static class Simulator
{
    public static void AssignMutationsRecursive(Clone currentClone, List<Clone> clones, AberrationsInfo aberrationsInfo,
        Random rnd, SimParams simParams)
    {
        foreach (int cloneId in currentClone.ChildrenIDs)
        {
            clones[cloneId].Karyotype = currentClone.SetKaryotype();
            for (int i = 0; i < clones[cloneId].MutCount; i++)
            {
                var aberration = aberrationsInfo.PickRandomMutation(rnd);
                string region = clones[cloneId].Karyotype.ApplyAberration(rnd, aberration,
                    aberrationsInfo.Map[aberration], clones[cloneId].Name);
                AssignFitness(clones[cloneId], simParams);
                new Abberation(clones[cloneId].Name, aberration.ToString(), region,
                    (float) clones[cloneId].deltaFitness);
            }

            AssignMutationsRecursive(clones[cloneId], clones, aberrationsInfo, rnd, simParams);
        }
    }

    public static List<Clone> GetClonePair(int distance, bool isFemale)
    {
        var parent = new Clone(0, -1, "1", 0, new Karyotype(isFemale));
        var child = new Clone(1, 0, "2", distance, new Karyotype(isFemale));
        parent.ChildrenIDs.Add(child.CloneId);
        return new List<Clone> {parent, child};
    }

    // Represents the limitation of space in the nucleus - more chromosomes ==> more stress
    // TODO: This needs to be validated
    private static float CalcStress(float stressFactor, int chromCount)
        => stressFactor * (float) Math.Pow(Math.Max(0, chromCount - 46), 2);

    public static void AssignFitness(Clone clone, SimParams simParams)
    {
        
        clone.deltaFitness = 0;
        float tsgOgFitness = 0;
        float stress = CalcStress(simParams.stressFraction, clone.Karyotype.ChromCount);

        var essentialStillThere = new List<Gen>();
        var tsgOgStillThere = new List<Gen>();
        foreach (var chr in clone.Karyotype.Chromosomes)
        {
            foreach (var region in chr._regions)
            {
                var chromNum = region.ChromId.ChromNum;
                essentialStillThere.AddRange(GenList.EssentialList[chromNum]
                    .FindAll(x => x.region.Start > region.Start && x.region.End < region.End));
                tsgOgStillThere.AddRange(GenList.TsgOgList[chromNum]
                    .FindAll(x => x.region.Start > region.Start && x.region.End < region.End));
            }
        }
        var essentialsLost = FindMissingGenes(essentialStillThere, GenList.EssentialList);
        var tsgOgLost = FindMissingGenes(tsgOgStillThere, GenList.TsgOgList)
                .Concat(tsgOgStillThere.GroupBy(x => x).Where(g => g.Count() < 2).Select(g => g.Key));
        var tsgOgAdded = tsgOgStillThere.GroupBy(x => x).Where(g => g.Count() > 2).Select(g => g.Key);
        float essentialityFitness = essentialsLost.Sum(essential => essential.deltaFitness);
        foreach (var tsgOg in tsgOgAdded)
        {
            tsgOgFitness -= tsgOg.deltaFitness;
        }
        foreach (var tsgOg in tsgOgLost)
        {
            tsgOgFitness += tsgOg.deltaFitness;
        }

        clone.deltaFitness = stress + (simParams.tsgOgFraction*tsgOgFitness) + 
            (simParams.essentialFraction*essentialityFitness);
    }

    private static List<Gen> FindMissingGenes(List<Gen> genes, Dictionary<ChromNum, List<Gen>> dict)
    {
        List<Gen> missingGenes = new List<Gen>();
        foreach (var list in dict)
        {
            missingGenes.AddRange(list.Value.Except(genes));
        }

        return missingGenes;
    }
}