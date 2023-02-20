using SimChA.DataTypes;
using ExtremeBinDist = Extreme.Statistics.Distributions.BinomialDistribution;

namespace SimChA.Simulation;

public class Simulator
{
    private readonly Random _rnd;
    private readonly SimParams _simParams;
    private readonly Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> _geneLists;
    
    public Simulator(
        Random rnd, 
        SimParams simParams, 
        Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> geneLists)
    {
        _rnd = rnd;
        _simParams = simParams;
        _geneLists = geneLists;
    }
    
    // calculate the number of nodes in the tree given by clones from the rootClone
    private static int GetTreeNodeCount(Clone root, List<Clone> clones)
        => 1 + root.ChildrenIDs.Select(id => GetTreeNodeCount(clones[id], clones)).Sum();

    public List<CNEvent> ApplyEvents(Clone rootClone, List<Clone> clones)
    {
        List<CNEvent> events = new();
        var progress = (1, GetTreeNodeCount(rootClone, clones) - 1);
        AssignMutationsRecursive(rootClone, clones, events, progress);
        return events;
    }
    
    private void AssignMutationsRecursive(Clone node, List<Clone> clones, List<CNEvent> events, (int, int) progress)
    {
        foreach (var child in node.ChildrenIDs.Select(cloneId => clones[cloneId]))
        {
            child.Karyotype = node.CopyKaryotype();
            double oldFitness = node.Karyotype.FitnessVal;
            int parentMutations = GetMutations(node, clones);
            for (int i = 0; i < child.DistToParent; i++)
            {
                Console.Write($"Clone {progress.Item1++}/{progress.Item2}, Mut {i+1}/{child.DistToParent}.\r");
                var sig = SignatureHelper.PickRandomSignature(_rnd, _simParams.Signatures);
                var eventP = SignatureHelper.PickRandomEventP(_rnd, sig.Events);
                string eventString = child.Karyotype.ApplyAberration(_rnd, eventP);
                double newFitness = child.Karyotype.UpdateFitness(_geneLists, _simParams.Fitness);
                int mutationCount = parentMutations + 1 + i;
                var abberation = new CNEvent(child.Name, eventP.Type, mutationCount,
                    eventString, newFitness - oldFitness, newFitness);
                events.Add(abberation);
                oldFitness = newFitness;
            }
            AssignMutationsRecursive(child, clones, events, progress);
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