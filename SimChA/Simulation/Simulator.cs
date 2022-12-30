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
    
    // calculate the number of nodes in the tree given by clones from the rootClone
    private static int GetTreeNodeCount(Clone root, List<Clone> clones)
        => 1 + root.ChildrenIDs.Select(id => GetTreeNodeCount(clones[id], clones)).Sum();

    public List<Abberation> AssignMutations(Clone rootClone, List<Clone> clones)
    {
        List<Abberation> abberationList = new();
        int numNodes = GetTreeNodeCount(rootClone, clones) - 1;
        int nodeNo = 1;
        AssignMutationsRecursive(rootClone, clones, abberationList, ref nodeNo, numNodes);
        return abberationList;
    }
    
    private void AssignMutationsRecursive(Clone node, List<Clone> clones, List<Abberation> abberationList, ref int cloneNo, int numNodes)
    {
        foreach (var child in node.ChildrenIDs.Select(cloneId => clones[cloneId]))
        {
            child.Karyotype = node.CopyKaryotype();
            float oldFitness = child.Fitness = node.Fitness;
            int parentMutations = GetMutations(node, clones);
            for (int i = 0; i < child.DistToParent; i++)
            {
                Console.Write($"Clone {cloneNo}/{numNodes}, Mut {i+1}/{child.DistToParent}.\r");
                var aberration = _aberrationsInfo.PickRandomMutation(_rnd);
                string eventString = child.Karyotype.ApplyAberration(_rnd, aberration, _aberrationsInfo.Map[aberration]);
                child.Fitness = CalcFitness(child.Karyotype);
                float deltaFitness = child.Fitness - oldFitness;
                oldFitness = child.Fitness;
                int mutationCount = parentMutations + 1 + i;
                var abberation = new Abberation(child.Name, aberration.ToString(), mutationCount,
                    eventString, deltaFitness, child.Fitness);
                abberationList.Add(abberation);
            }
            AssignMutationsRecursive(child, clones, abberationList, ref cloneNo, numNodes);
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
        var parent = new Clone(0, -1, "1", 0, new Karyotype(isFemale), 0);
        var child = new Clone(1, 0, "2", distance, new Karyotype(isFemale), distance);
        parent.ChildrenIDs.Add(child.CloneId);
        return new List<Clone> { parent, child };
    }

    // Represents the limitation of space in the nucleus - more chromosomes ==> more stress
    // TODO: This needs to be validated
    private static float CalcStress(float stressFactor, int chromCount)
        => stressFactor * (float)Math.Pow(Math.Max(0, chromCount - 46), 2);

    private float CalcFitness(Karyotype kar)
    {
        float stress = CalcStress(_simParams.StressFraction, kar.ChromCount);
        var essentialFound = kar.GetPresentGenes(_essentialGenes);
        var tsgOgFound = kar.GetPresentGenes(_tsgOgGenes);
        var essentialsMissing = FindMissingGenes(essentialFound, _essentialGenes);
        var tsgOgMissing = FindMissingGenes(tsgOgFound, _tsgOgGenes);
        var tsgOgCounts = tsgOgFound.GroupBy(x => x).ToList();
        float essentialityFitness = essentialsMissing.Sum(g => g.DeltaFitness);
        // Twice the value for missing genes (ploidy 0), -1 multiplicative factor for each missing gene (ploidy 1),
        // n - 2 for each overrepresented gene (ploidy 2+)
        float tsgOgFitness = 2*tsgOgMissing.Sum(g => g.DeltaFitness) 
                             -tsgOgCounts.Sum(g => g.Key.DeltaFitness * (g.Count() - 2));
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