using SimChA.DataTypes;
using EDists = Extreme.Statistics.Distributions;
using SimChA.Computation;

namespace SimChA.Simulation;

public class Simulator
{
    private readonly Random _rnd;
    private readonly SimParams _simParams;
    private readonly Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> _geneLists;

    private Signature _sig;
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
        _sig = SignatureHelper.PickRandomSignature(_rnd, _simParams.Signatures);
        // Return the event weight
        return SignatureHelper.PickRandomEventP(_rnd, _sig.Events);
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

    // The conditional probability of a given event occuring, given the signature
    private double MutationPotential(List<CNEventP> eventPs)
    {
        // Probability of picking this signature
        var signatures = _simParams.Signatures;
        var sig_potential = 1.0;
        if (signatures != null)
        {
            sig_potential = _sig.Prob / signatures.Sum(s => s.Prob);
        }
        // Probability of picking this set of events
        double event_potential_total = 1.0;
        foreach (CNEventP eP in eventPs)
        {
            event_potential_total *= eP.Prob;
        }
        // Normalize the event potential
        event_potential_total /= Math.Pow(_sig.Events.Sum(e => e.Prob), eventPs.Count);
        return sig_potential * event_potential_total;
    }

    private double Potential(Clone node, List<CNEventP> eventPs)
    {
        double fitness_potential = 0.0;

        // Probability of picking this set of events
        double event_potential_total = 1.0;
        double meanFitness = 0.0;

        // Make a temporary copy of the adult clone
        var karyotype = node.CopyKaryotype();
        // Apply the events
        foreach (var eventP in eventPs)
        {
            string eventString = karyotype.ApplyCNEvent(_rnd, eventP);
            // Update the probability for the event potential
            event_potential_total *= eventP.Prob;
        }

        double newFitness = karyotype.UpdateFitness(_geneLists, _simParams.Fitness);
        // Normalize the event potential
        event_potential_total /= Math.Pow(_sig.Events.Sum(e => e.Prob), eventPs.Count);

        double thetaFitness = 1.0;
        // Fitness potential is an exponential - exp[-theta * |fit - mean_fit|]
        fitness_potential = Math.Pow(newFitness, thetaFitness);
        //double mean_fit_potential = Exp();
        
        
        return fitness_potential * event_potential_total;
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
            float SwapEventP = 1.0f;
            float AlterEventStart = 0.5f;
            float AlterEventLength = 0.5f;

            Console.WriteLine("Generating starting events:");
            // Generate a starting set of mutations
            List<CNEventP> currentEventPs  = InitEvents(child.DistToParent);

            // Calculate the overall fitness of this clone
            double currentPotential = 0.0f;
            
            //string eventString = child.Karyotype.ApplyCNEvent(_rnd, eventP);
            //double newFitness = Fitness.Calculate(this, geneLists, fParams);

            // Now we perform the Metropolis-Hastings algorithm
            Console.WriteLine("\nPerforming Metropolis-Hastings");
            for (int i = 0; i < n_samples; i++)
            {   
                var proposedEventPs = currentEventPs;
                //Console.WriteLine($"{i%1000}");
                // Printing out benchmarks in sampling
                if (i > burn_in && i%1000 == 0)
                {
                    Console.Write($"\rEvent {i-burn_in}   ");
                }
                
                // Select a random CNEventP to modify
                int index = _rnd.Next(proposedEventPs.Count);
                // Choose whether to swap the event entirely
                if (EDists.ContinuousUniformDistribution.Sample(_rnd, 0, 1) < SwapEventP)
                {
                    proposedEventPs[index] = newEventP();
                }
                // With the newly selected event, we need to calculate the new
                // fitness of the clone
                // TODO: calc fitness
                double proposalPotential = MutationPotential(proposedEventPs);//EDists.ContinuousUniformDistribution.Sample(_rnd, 0, 1);
                double acceptProb = proposalPotential/currentPotential;
                if (acceptProb >= 1 || acceptProb > EDists.ContinuousUniformDistribution.Sample(_rnd,0,1))
                {
                    currentPotential = proposalPotential;
                    currentEventPs   = proposedEventPs;
                }

            }
            //counter++;
            //ApplyCNEventsRec(child, clones, events, ref counter);
        }
    }
}