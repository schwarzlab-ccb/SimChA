using SimChA.DataTypes;
using EDists = Extreme.Statistics.Distributions;
using SimChA.Computation;

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
            for (int mutNo = 0; mutNo < child.DistToParent; mutNo++)
            {
                Console.Write($"\rClone {counter}/{clones.Count-1}. Event {mutNo+1}/{child.DistToParent}.");
                var sig = SignatureHelper.PickRandomSignature(_rnd, _simParams.Signatures);
                var eventP = SignatureHelper.PickRandomEventP(_rnd, sig.Events);
                string eventString = child.Karyotype.ApplyCNEvent(_rnd, eventP);
                double newFitness = child.Karyotype.UpdateFitness(_geneLists, _simParams.Fitness);
                int mutationCount = parentMutations + 1 + mutNo;
                double dFit = newFitness - oldFitness;
                var abberation = new CNEvent(child.CloneId, mutationCount, eventP.Type, eventString, dFit, newFitness);
                events.Add(abberation);
                oldFitness = newFitness;
            }
            counter++;
            ApplyCNEventsRec(child, clones, events, ref counter);
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

    public void MCSampleEvents(Clone rootClone, List<Clone> clones, float fitness)
    {
        List<CNEvent> events = new();
        int counter = 1;
        MCSampleCNEventsRec(rootClone, clones, events, fitness, ref counter);
        Console.WriteLine();
        return;
    }


    private CNEventP newEventP()
    {
        // Pick a random allowed signature
        var sig = SignatureHelper.PickRandomSignature(_rnd, _simParams.Signatures);
        // Return the event weight
        return SignatureHelper.PickRandomEventP(_rnd, sig.Events);
    }

    // Method to generate the starting events
    private List<CNEventP> InitEvents(int nMutations)
    {
        List<CNEventP> tempEventPs = new();
        for (int mutNo = 0; mutNo < nMutations; mutNo++)
        {   
            tempEventPs.Add(newEventP());
        }
        return tempEventPs;
    }

    private void MCSampleCNEventsRec(Clone node, List<Clone> clones, List<CNEvent> events, float fitness, ref int counter)
    {
        foreach (var child in node.ChildrenIDs.Select(cloneId => clones[cloneId]))
        {
            child.Karyotype = node.CopyKaryotype();
            double oldFitness = node.Karyotype.FitnessVal;
            int parentMutations = GetMutations(node, clones);
            // TODO: Maybe I need a dummy Karyotype object to apply events to

            // Parameters needed for the MH algorithm
            // Number of trial events
            int burn_in = 1000;
            int n_samples = burn_in + 20000;
            float SwapEventP = 0.5f;
            float AlterEventStart = 0.5f;
            float AlterEventLength = 0.5f;

            Console.WriteLine("Generating starting events:");
            // Generate a starting set of mutations
            List<CNEventP> nowEventPs  = InitEvents(child.DistToParent);
            // Storage for the set of events giving best agreement with mean fitness
            List<CNEventP> bestEventPs = nowEventPs;
            
            // Calculate the overall fitness of this clone
            double nowFitness = 0.0f;
            double bestFitness  = nowFitness;
            //string eventString = child.Karyotype.ApplyCNEvent(_rnd, eventP);
            //double newFitness = Fitness.Calculate(this, geneLists, fParams);

            // Now we perform the Metropolis-Hastings algorithm
            Console.WriteLine("\nPerforming Metropolis-Hastings");
            for (int i = 0; i < n_samples; i++)
            {   
                //Console.WriteLine($"{i%1000}");
                // Printing out benchmarks in sampling
                if (i > burn_in && i%1000 == 0)
                {
                    Console.Write($"\rEvent {i-burn_in}   ");
                }
                // Select a random CNEventP to modify
                int index = _rnd.Next(bestEventPs.Count);
                // Choose whether to swap the event entirely
                if (EDists.ContinuousUniformDistribution.Sample(_rnd, 0, 1) < SwapEventP)
                {
                    nowEventPs[index] = newEventP();
                }


            }
            //counter++;
            //ApplyCNEventsRec(child, clones, events, ref counter);
        }
    }
}