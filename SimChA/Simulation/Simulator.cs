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
        int counter = 1;
        ApplyCNEventsRec(rootClone, clones, events, ref counter);
        Console.WriteLine();
        return events;
    }
    
    private void ApplyCNEventsRec(Clone node, List<Clone> clones, List<CNEvent> events, ref int counter)
    {
        foreach (var child in node.ChildrenIDs.Select(cloneId => clones[cloneId]))
        {
            child.Karyotype = node.CopyKaryotype();
            double oldFitness = node.Karyotype.FitnessVal;
            int parentMutations = GetMutations(node, clones);
            int i = 0;
            int j = 0;
            while (i < child.DistToParent)
            {
                Console.Write($"\rClone {counter}/{clones.Count-1}. Event {i+1}/{child.DistToParent}. Attempt {++j}");
                var sig = SignatureHelper.PickRandomSignature(_rnd, _simParams.Signatures);
                var eventP = SignatureHelper.PickRandomEventP(_rnd, sig.Events);
                string eventString = child.Karyotype.ApplyAberration(_rnd, eventP);
                double newFitness = child.Karyotype.UpdateFitness(_geneLists, _simParams.Fitness);
                int mutationCount = parentMutations + 1 + i;
                double dFit = newFitness - oldFitness;
                var abberation = new CNEvent(child.Name, eventP.Type, mutationCount, eventString, dFit, newFitness);
                events.Add(abberation);
                oldFitness = newFitness;
                i++;
                j = 0;
            }
            counter++;
            ApplyCNEventsRec(child, clones, events, ref counter);
        }
    }
    
    private static int GetMutations(Clone clone, List<Clone> clones)
    {
        int mutCount = clone.ParentId != -1 ? GetMutations(clones[clone.ParentId], clones) : 0;
        mutCount += clone.DistToParent;
        return mutCount;
    }

    public static List<Clone> MakeClones(int distance, int repeats, bool sexXX)
    {
        var parent = new Clone(0, -1, "0-0", 0, new Karyotype(sexXX), 0);
        var clones = new List<Clone> { parent };
        for (var i = 1; i <= repeats; i++)
        {
            var child = new Clone(i, 0, $"{i}-{distance}", distance, new Karyotype(sexXX), distance);
            parent.ChildrenIDs.Add(child.CloneId);
            clones.Add(child);
        }
        return clones;
    }

    public static List<Clone> ClonesFromProfiles(Dictionary<string, Karyotype> profiles)
    {
        var i = 1;
        var res = profiles.Select(pair => new Clone(i++, 0, pair.Key, 1, pair.Value, 1));
        return res.ToList();
    }
}