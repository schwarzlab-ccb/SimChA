using SimChA.DataTypes;
using ExtremeBinDist = Extreme.Statistics.Distributions.BinomialDistribution;
using SimChA.IO;

namespace SimChA.Simulation;

public class Simulator
{
    private readonly AberrationsInfo _aberrationsInfo;
    private readonly Random _rnd;
    private readonly SimParams _simParams;
    private readonly Dictionary<ChrNo, List<Gene>> _tsgOgGenes;
    private readonly Dictionary<ChrNo, List<Gene>> _essentialGenes;
    
    public Simulator(
        AberrationsInfo aberrationsInfo, 
        Random rnd, 
        SimParams simParams, 
        Dictionary<ChrNo, List<Gene>> tsgOgGenes, 
        Dictionary<ChrNo, List<Gene>> essentialGenes)
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
            double oldFitness = node.Karyotype.FitnessVal;
            int parentMutations = GetMutations(node, clones);
            for (int i = 0; i < child.DistToParent; i++)
            {
                Console.Write($"Clone {cloneNo}/{numNodes}, Mut {i+1}/{child.DistToParent}.\r");
                var aberration = _aberrationsInfo.PickRandomMutation(_rnd);
                string eventString = child.Karyotype.ApplyAberration(_rnd, aberration, _aberrationsInfo.Map[aberration]);
                double newFitness = child.Karyotype.UpdateFitness(_essentialGenes, _tsgOgGenes, _simParams);
                int mutationCount = parentMutations + 1 + i;
                var abberation = new Abberation(child.Name, aberration, mutationCount,
                    eventString, newFitness - oldFitness, newFitness);
                abberationList.Add(abberation);
                oldFitness = newFitness;
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
}