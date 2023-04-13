using SimChA.DataTypes;

namespace SimChA.Simulation;

public class Simulator
{
    private readonly Random _rnd;
    private readonly SimParams _simParams;
    private readonly Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> _geneLists;

    private MCParams? McParams => _simParams.MCParams;
    private List<Signature>? Signatures => _simParams.Signatures;
    
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
                var sig = SignatureHelper.RndSignature(_rnd, Signatures);
                var eventP = SignatureHelper.RndEventP(_rnd, sig.Events);
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

    public static List<Clone> MakeClones(Random rnd, int repeats, bool sexXX, int distance, Distribution distribution)
    {
        var parent = new Clone(0, -1, "0-0", 0, new Karyotype(sexXX), 0);
        var clones = new List<Clone> { parent };
        for (var i = 1; i <= repeats; i++)
        {
            double sample = Sampling.SampleDist(rnd, distribution);
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

    public List<CNEvent> MCSampleEvents(Clone rootClone, List<Clone> clones, Dictionary<string, double> fitnessMap)
    {
        List<CNEvent> events = new();
        var counter = 1;
        Console.WriteLine("\nPerforming Metropolis-Hastings");
        MCSampleCNEventsRec(rootClone, clones, fitnessMap, events, ref counter);
        Console.WriteLine();
        return events;
    }

    public List<CNEventProperties> InitEvents(Clone node, Signature sig, int nMutations)
    {
        var eventPs = InitEventPs(sig, nMutations).ToList();
        List<CNEventProperties> eventProperties = new List<CNEventProperties>();
        var karyotype = node.CopyKaryotype();
        foreach (var e in eventPs)
            eventProperties.Add(karyotype.GenerateCNEventProperties(_rnd, e));

        return eventProperties;
    }
    public IEnumerable<CNEventP> InitEventPs(Signature sig, int nMutations)
        => Enumerable.Range(0, nMutations).Select(_ => SignatureHelper.RndEventP(_rnd, sig.Events));

    // The conditional probability of this set of events occuring, 
    // given the individual events and the signature
    // TODO: do we need to change the eventPs as a result
    private double Potential(Clone node, Signature sig, List<CNEventProperties> events)
    {
        // Probability of picking this set of events
        double eventPotentialTotal = 1.0;
        double meanFitness = 0.0; //_fitnessDict[node.CloneId.ToString()];
        // Probability of picking this signature, no need to normalize, since it will
        // divide out when comparing sets of events
        var sigPotential = sig.Prob;
        
        /* Currently assume one signature
        if (signatures != null)
        {
            sigPotential = sig.Prob / signatures.Sum(s => s.Prob);
        }*/
        
        var karyotype = node.CopyKaryotype();
        foreach (var e in events)
        {
            karyotype.ApplyCNEvent(_rnd, e.EventP);
            eventPotentialTotal *= e.EventP.Prob;
        }

        double newFitness = karyotype.UpdateFitness(_geneLists, _simParams.Fitness);
        // Normalize the event potential
        //eventPotentialTotal /= Math.Pow(sig.Events.Sum(e => e.Prob), eventPs.Count);
        double dFit = newFitness - meanFitness;
        // Fitness potential is an exponential - exp[-theta * |fit - mean_fit|]
        double fitnessPotential = Math.Exp(-McParams.ThetaFitness * Math.Abs(dFit));
        // Gaussian form
        //double fitnessPotential = Math.Exp(-McParams.ThetaFitness*Math.Pow(dFit,2));
        
        return fitnessPotential * eventPotentialTotal * sigPotential;
    }

    private void MCSampleCNEventsRec(Clone node, List<Clone> clones, Dictionary<string, double> fitnessMap, List<CNEvent> events,  ref int counter)
    {
        foreach (var child in node.ChildrenIDs.Select(cloneId => clones[cloneId]))
        {
            // Skip any clones that don't have any mutational distance from their
            // parent
            if (child.DistToParent == 0)
                continue;
            
            // Initialize all the relevant quantities
            double oldFitness = node.Karyotype.FitnessVal;
            int parentMutations = GetMutations(node, clones);
            
            // TODO: Maybe I need a dummy Karyotype object to apply events to
            // Parameters needed for the MH algorithm
            // Number of trial events
            //int nSamples = _mcParams.NumBurnIn + _mcParams.NumSamples;
            float alterEventStart = 0.5f;
            float alterEventLength = 0.5f;

            // Generate a starting set of mutations and its potential
            var sig = SignatureHelper.RndSignature(_rnd, Signatures);
            var currentEventProps = InitEvents(node, sig, child.DistToParent);
            double currentPotential = Potential(node, sig, currentEventProps);

            var startEventProps = currentEventProps.ToList();
            // Now we perform the Metropolis-Hastings algorithm
            // and sample a set of events that give the closest agreement with
            // fitness given by SMITH
            for (int i = 0; i < McParams.NumSamples; i++)
            {
                var proposedEventProps = currentEventProps.ToList();
                // Select a random CNEventP to modify
                int index = _rnd.Next(proposedEventProps.Count);
                // Choose whether to swap the event entirely
                if (_rnd.NextDouble() < McParams.SwapEventP)
                {
                    var e = SignatureHelper.RndEventP(_rnd, sig.Events);
                    proposedEventProps[index] = node.Karyotype.GenerateCNEventProperties(_rnd, e);
                }
                // Otherwise we modify some quantity of the event
                // TODO: Implement
                else
                {
                    var eventProp = proposedEventProps[index];
                }
                // With the newly selected event, we need to calculate the new
                // fitness of the clone
                double proposalPotential = Potential(node, sig, currentEventProps);
                double acceptProb = proposalPotential/currentPotential;
                
                if (acceptProb < 1 && acceptProb <= _rnd.NextDouble()) 
                    continue;
                
                currentPotential = proposalPotential;
                currentEventProps = proposedEventProps;
            }

            // Finalize the mutated karyotype by applying the best-fit set of events
            // and move onto the next clone
            // Copy its parent
            child.Karyotype = node.CopyKaryotype();
            for (int mutNo = 0; mutNo < currentEventProps.Count; mutNo++)
            {
                Console.Write($"\rClone {counter}/{clones.Count-1}. Event {mutNo+1}/{child.DistToParent}.");
                var eventProperties = currentEventProps[mutNo];
                string eventString = child.Karyotype.ApplyEventProperties(_rnd, eventProperties);
                double newFitness = child.Karyotype.UpdateFitness(_geneLists, _simParams.Fitness);
                int mutationCount = parentMutations + 1 + mutNo;
                double dFit = newFitness - oldFitness;
                var abberation = new CNEvent(child.CloneId, mutationCount, eventProperties.EventType, eventString, dFit, newFitness);
                events.Add(abberation);
            }
            
            counter++;
            MCSampleCNEventsRec(child, clones, fitnessMap, events, ref counter);
        }
    }
}