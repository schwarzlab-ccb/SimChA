using SimChA.DataTypes;
using ExtremeBinDist = Extreme.Statistics.Distributions.BinomialDistribution;
using SimChA.IO;

namespace SimChA.Simulation;

public class Simulator
{
    private readonly AberrationsInfo _aberrationsInfo;
    private readonly Random _rnd;
    private readonly SimParams _simParams;
    private readonly Dictionary<ChromNum, List<Gene>> _tsgOgGenes;
    private readonly Dictionary<ChromNum, List<Gene>> _essentialGenes;
    
    public Simulator(
        AberrationsInfo aberrationsInfo, 
        Random rnd, 
        SimParams simParams, 
        Dictionary<ChromNum, List<Gene>> tsgOgGenes, 
        Dictionary<ChromNum, List<Gene>> essentialGenes)
    {
        _aberrationsInfo = aberrationsInfo;
        _rnd = rnd;
        _simParams = simParams;
        _tsgOgGenes = tsgOgGenes;
        _essentialGenes = essentialGenes;
    }
    
    public List<Abberation> AssignMutations(Clone rootClone, List<Clone> clones)
    {
        List<Abberation> abberationList = new();
        AssignMutationsRecursive(rootClone, clones, abberationList);
        return abberationList;
    }
    
    private void AssignMutationsRecursive(Clone node, List<Clone> clones, List<Abberation> abberationList)
    {
        foreach (var child in node.ChildrenIDs.Select(cloneId => clones[cloneId]))
        {
            child.Karyotype = node.CopyKaryotype();
            child.Fitness = node.Fitness;
            int parentMutations = GetMutations(node, clones);
            for (int i = 0; i < child.DistToParent; i++)
            {
                Console.Write($"Clone {child.CloneId}, Mut {i+1}/{child.DistToParent}.\r");
                float oldFitness = child.Fitness;
                var aberration = _aberrationsInfo.PickRandomMutation(_rnd);
                string eventString = child.Karyotype.ApplyAberration(_rnd, aberration, _aberrationsInfo.Map[aberration]);
                child.Fitness = CalcFitness(child);
                int mutationCount = parentMutations + 1 + i;
                float deltaFitness = child.Fitness - oldFitness;
                var abberation = new Abberation(child.Name, aberration.ToString(), mutationCount,
                    eventString, deltaFitness, child.Fitness);
                abberationList.Add(abberation);
            }
            AssignMutationsRecursive(child, clones, abberationList);
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

    private float CalcFitness(Clone clone)
    {
        float stress = CalcStress(_simParams.StressFraction, clone.Karyotype.ChromCount);
        List<Gene> essentialFound = new();
        List<Gene> tsgOgFound = new();
        foreach (var chr in clone.Karyotype.Chromosomes)
        {
            foreach (var region in chr.GetAllRegions())
            {
                var chromNum = region.ChromId.ChromNum;
                // TODO use overlap
                essentialFound.AddRange(_essentialGenes[chromNum]
                    .FindAll(x => x.Region.Start > region.Start && x.Region.End < region.End));
                tsgOgFound.AddRange(_tsgOgGenes[chromNum]
                    .FindAll(x => x.Region.Start > region.Start && x.Region.End < region.End));
            }
        }
        var essentialsMissing = FindMissingGenes(essentialFound, _essentialGenes);
        var tsgOgMissing = FindMissingGenes(tsgOgFound, _tsgOgGenes);
        var tsgOgCounts = tsgOgFound.GroupBy(x => x).ToList();
        var tsgOgInsufficient = tsgOgCounts.Where(g => g.Count() < 2).Select(g => g.Key);
        // TODO: This is probably not correct! The delta will be the same for ploidy 0 and 1!
        var tsgOgLost = tsgOgMissing.Concat(tsgOgInsufficient);
        var tsgOgDuplicated = tsgOgCounts.Where(g => g.Count() > 2).Select(g => g.Key);
        float essentialityFitness = essentialsMissing.Sum(g => g.DeltaFitness);
        float tsgOgFitness = tsgOgLost.Sum(g => g.DeltaFitness) - tsgOgDuplicated.Sum(g => g.DeltaFitness);

        // parametrized linear combination of factors
        return stress + _simParams.TsgOgFraction * tsgOgFitness + _simParams.EssentialFraction * essentialityFitness;
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