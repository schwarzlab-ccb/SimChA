using SimChA.DataTypes;
using EDists = Extreme.Statistics.Distributions;

namespace SimChA.Simulation;

public class Simulator
{
    private readonly Random _rnd;
    private readonly FitnessParams _fitness;
    private readonly List<Signature> _signatures;
    private readonly Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> _geneLists;
    
    
    public Simulator(
        Random rnd, 
        FitnessParams fitnessParams,
        List<Signature> signatures,
        Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> geneLists)
    {
        _rnd = rnd;
        _fitness = fitnessParams;
        _signatures = signatures;
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
    
    private void ApplyCNEventsRec(Clone node, List<Clone> clones, List<CNEvent> eventSeq, ref int counter)
    {
        foreach (var child in node.ChildrenIDs.Select(cloneId => clones[cloneId]))
        {
            child.Karyotype = node.CopyKaryotype();
            double oldFitness = node.Karyotype.FitnessVal;
            int parentMutations = GetMutations(node, clones);
            for (var mutNo = 0; mutNo < child.DistToParent; mutNo++)
            {
                int mixtureIndex = Sampling.PickRandomIndex(_rnd, child.SigMixture);
                var cnEventParams = _signatures[mixtureIndex].Events;
                Console.Write($"\rClone {counter}/{clones.Count-1}. Event {mutNo+1}/{child.DistToParent}.");
                var eventP = Sampling.PickRandomEventP(_rnd, cnEventParams);
                string eventString = child.Karyotype.ApplyCNEvent(_rnd, eventP);
                double newFitness = child.Karyotype.UpdateFitness(_geneLists, _fitness);
                int mutationCount = parentMutations + 1 + mutNo;
                double dFit = newFitness - oldFitness;
                var abberation = new CNEvent(child.CloneId, mutationCount, eventP.Type, eventString, dFit, newFitness);
                eventSeq.Add(abberation);
                oldFitness = newFitness;
            }
            counter++;
            ApplyCNEventsRec(child, clones, eventSeq, ref counter);
        }
    }
    
    private static int GetMutations(Clone clone, IReadOnlyList<Clone> clones)
    {
        int mutCount = clone.ParentId != -1 ? GetMutations(clones[clone.ParentId], clones) : 0;
        mutCount += clone.DistToParent;
        return mutCount;
    }
    
    private static double GetSample(Random rnd, Distribution dist)
    {
        return dist switch
        {
            Distribution.Exponential => EDists.ExponentialDistribution.Sample(rnd, 1),
            Distribution.Normal => EDists.NormalDistribution.Sample(rnd, 1, 1),
            _ => 1
        };
    }

    public static List<Clone> MakeClones(Random rnd, int repeats, bool sexXX, int distance, Distribution distribution)
    {
        var parent = new Clone(0, -1, "0-0", 0, new Karyotype(sexXX), 0);
        var clones = new List<Clone> { parent };
        for (var i = 1; i <= repeats; i++)
        {
            double sample = GetSample(rnd, distribution);
            var mutCount = (int) Math.Round(distance * sample);
            var child = new Clone(i, 0, $"{i}-{mutCount}", mutCount, new Karyotype(sexXX), mutCount);
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