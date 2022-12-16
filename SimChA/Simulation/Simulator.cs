using Extreme.Statistics.Distributions;
using SimChA.DataTypes;
using ExtremeBinDist = Extreme.Statistics.Distributions.BinomialDistribution;
using SimChA.IO;

namespace SimChA.Simulation;

public static class Simulator
{
    public static void AssignMutationsRecursive(Clone currentClone, List<Clone> clones, List<Abberation> abberationList,
        AberrationsInfo aberrationsInfo, Random rnd, SimParams simParams)
    {
        foreach (var child in currentClone.ChildrenIDs.Select(cloneId => clones[cloneId]))
        {
            child.Karyotype = currentClone.CopyKaryotype();
            child.Fitness = currentClone.Fitness;
            int parentMutations = GetMutations(clones[currentClone.CloneId], clones);
            for (int i = 0; i < child.DistToParent; i++)
            {
                float oldFitness = child.Fitness;
                var aberration = aberrationsInfo.PickRandomMutation(rnd);
                string eventString = child.Karyotype.ApplyAberration(rnd, aberration, aberrationsInfo.Map[aberration]);
                child.Fitness = CalcFitness(child, simParams);
                int mutationCount = parentMutations + 1 + i;
                float deltaFitness = child.Fitness - oldFitness;
                var abberation = new Abberation(child.Name, aberration.ToString(), mutationCount,
                    eventString, deltaFitness, child.Fitness);
                abberationList.Add(abberation);
            }
            AssignMutationsRecursive(child, clones, abberationList, aberrationsInfo, rnd, simParams);
        }
    }

    private static int GetMutations(Clone clone, List<Clone> clones)
    {
        int mutCount = clone.ParentId != -1 ? GetMutations(clones[clone.ParentId], clones) : 0;
        mutCount += clone.DistToParent;
        return mutCount;
    }

    public static List<Clone> MakeClonePair(int distance, bool isFemale)
    {
        var parent = new Clone(0, -1, "1", 0, new Karyotype(isFemale));
        var child = new Clone(1, 0, "2", distance, new Karyotype(isFemale));
        parent.ChildrenIDs.Add(child.CloneId);
        return new List<Clone> { parent, child };
    }

    // Represents the limitation of space in the nucleus - more chromosomes ==> more stress
    // TODO: This needs to be validated
    private static float CalcStress(float stressFactor, int chromCount)
        => stressFactor * (float)Math.Pow(Math.Max(0, chromCount - 46), 2);

    private static float CalcFitness(Clone clone, SimParams simParams)
    {
        float stress = CalcStress(simParams.StressFraction, clone.Karyotype.ChromCount);
        List<Gene> essentialFound = new();
        List<Gene> tsgOgFound = new();
        foreach (var chr in clone.Karyotype.Chromosomes)
        {
            foreach (var region in chr.GetAllRegions())
            {
                var chromNum = region.ChromId.ChromNum;
                // TODO use overlap
                essentialFound.AddRange(GeneList.EssentialList[chromNum]
                    .FindAll(x => x.Region.Start > region.Start && x.Region.End < region.End));
                tsgOgFound.AddRange(GeneList.TsgOgList[chromNum]
                    .FindAll(x => x.Region.Start > region.Start && x.Region.End < region.End));
            }
        }
        var essentialsMissing = FindMissingGenes(essentialFound, GeneList.EssentialList);
        var tsgOgMissing = FindMissingGenes(tsgOgFound, GeneList.TsgOgList);
        var tsgOgCounts = tsgOgFound.GroupBy(x => x).ToList();
        var tsgOgInsufficient = tsgOgCounts.Where(g => g.Count() < 2).Select(g => g.Key);
        // TODO: This is probably not correct! The delta will be the same for ploidy 0 and 1!
        var tsgOgLost = tsgOgMissing.Concat(tsgOgInsufficient);
        var tsgOgDuplicated = tsgOgCounts.Where(g => g.Count() > 2).Select(g => g.Key);
        float essentialityFitness = essentialsMissing.Sum(g => g.DeltaFitness);
        float tsgOgFitness = tsgOgLost.Sum(g => g.DeltaFitness) - tsgOgDuplicated.Sum(g => g.DeltaFitness);

        // parametrized linear combination of factors
        return stress + simParams.TsgOgFraction * tsgOgFitness + simParams.EssentialFraction * essentialityFitness;
    }

    private static List<Gene> FindMissingGenes(List<Gene> presentGenes, Dictionary<ChromNum, List<Gene>> geneList)
    {
        var missingGenes = new List<Gene>();
        foreach (var (_, allGenes) in geneList)
        {
            missingGenes.AddRange(allGenes.Except(presentGenes));
        }
        return missingGenes;
    }
}