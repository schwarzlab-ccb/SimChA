using SimChA.DataTypes;
using SimChA.EventData;
using EDists = Extreme.Statistics.Distributions;

namespace SimChA.Simulation;

public class Simulator
{
    private readonly Random _rnd;
    private readonly FitnessParams _fitness;
    private readonly List<Signature> _signatures;
    private readonly MCParams _mcParams;
    private readonly Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> _geneLists;

    public List<Signature> SelectedSignatures;

    public Simulator(
        Random rnd,
        FitnessParams fitnessParams,
        List<Signature> signatures,
        MCParams mcParams,
        Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> geneLists)
    {
        _rnd = rnd;
        _fitness = fitnessParams;
        _signatures = signatures;
        _mcParams = mcParams;
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
                Console.Write($"\rClone {counter}/{clones.Count - 1}. Event {mutNo + 1}/{child.DistToParent}.");
                var eventP = Sampling.PickRandomEventP(_rnd, cnEventParams);
                var eventData = Sampling.GenerateCNEventData(_rnd, child.Karyotype, eventP);
                eventData.ApplyEvent(child.Karyotype);
                double newFitness = child.Karyotype.UpdateFitness(_geneLists, _fitness);
                int mutationCount = parentMutations + 1 + mutNo;
                double dFit = newFitness - oldFitness;
                var abberation = new CNEvent(child.CloneId, mutationCount, eventP.Type, eventData.ToString(), dFit, newFitness);
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
    
    public static List<Clone> MakeClones(Random rnd, int repeats, bool sexXX, int distance, Distribution distribution)
    {
        var parent = new Clone(0, -1, "0-0", 0, new Karyotype(sexXX), 0);
        var clones = new List<Clone> { parent };
        for (var i = 1; i <= repeats; i++)
        {
            double sample = Sampling.SampleDist(rnd, distribution);
            var mutCount = (int)Math.Round(distance * sample);
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
        MCSampleCNEventsRec(rootClone, clones, fitnessMap, events, ref counter);
        return events;
    }

    public List<BaseEventData> InitEvents(Clone node, int nMutations)
    {
        var eventPs = InitEventPs(node, nMutations);
        return eventPs.Select(e => Sampling.GenerateCNEventData(_rnd, node.Karyotype, e)).ToList();
    }

    public List<CNEventP> InitEventPs(Clone node, int nMutations)
    {
        List<CNEventP> eventPs = new List<CNEventP>();
        // Reset the selected signatures
        SelectedSignatures = new List<Signature>();
        for (int i = 0; i < nMutations; i++)
        {
            int mixtureIndex = Sampling.PickRandomIndex(_rnd, node.SigMixture);
            var sig = _signatures[mixtureIndex];
            SelectedSignatures.Add(sig);
            var eventP = Sampling.PickRandomEventP(_rnd, sig.Events);
            eventPs.Add(eventP);
        }
        return eventPs;
    }

    // The conditional probability of this set of events occuring, 
    // given the individual events and the signature
    public (double,double) LogPotential(Clone node, Dictionary<string, double> fitnessMap,
        List<BaseEventData> events, ref bool thresholdAccept)
    {
        double eventPotentialTotal = 0.0;
        double targetFitness = fitnessMap[node.Name];

        // Create a dummy karyotype for the events to act on
        var karyotype = node.CopyKaryotype();
        double sigPotential = 0.0;

        // Probability of picking each event and their corresponding signature
        for (int i = 0; i < events.Count; i++)
        {
            sigPotential += Math.Log(SelectedSignatures[i].Prob) - Math.Log(_signatures.Sum(sig => sig.Prob));
            events[i].ApplyEvent(karyotype);
            eventPotentialTotal += Math.Log(events[i].EventP.Prob) - Math.Log(SelectedSignatures[i].Events.Sum(e => e.Prob));
        }
        
        double newFitness = karyotype.UpdateFitness(_geneLists, _fitness);
        double dFit = newFitness - targetFitness;
        thresholdAccept = Math.Abs(dFit / targetFitness) < _mcParams.ThresholdFit;

        // Fitness potential is an exponential - exp[-theta * |fit - mean_fit|]
        double fitnessPotential = Math.Exp(-_mcParams.ThetaFitness * Math.Abs(dFit));
        // Gaussian form
        //double fitnessPotential = Math.Exp(-_mcParams.ThetaFitness*Math.Pow(dFit,2));

        return (Math.Log(fitnessPotential) + eventPotentialTotal + sigPotential, newFitness);
    }
    private void MCSampleCNEventsRec(Clone node, List<Clone> clones, Dictionary<string, double> fitnessMap,
        List<CNEvent> events, ref int counter)
    {
        foreach (var child in node.ChildrenIDs.Select(cloneId => clones[cloneId]))
        {
            // Skip
            Console.WriteLine($"\nTarget Fitness: {fitnessMap[child.Name]}");
            // Initialize all the relevant quantities
            double oldFitness = node.Karyotype.FitnessVal;
            int parentMutations = GetMutations(node, clones);
            child.Karyotype = node.CopyKaryotype();
            
            // Skip children with no mutational distance
            if (child.DistToParent == 0)
            {
                counter++;
                MCSampleCNEventsRec(child, clones, fitnessMap, events, ref counter);
                continue;
            }
            // Parameters needed for the MH algorithm
            float alterEventStart = 0.5f;
            float alterEventLength = 0.5f;
            // The tracker to accept sets of events if they are within the threshold tolerance from
            // the target fitness
            bool thresholdAccept = false;

            // Generate a starting set of mutations and its potential
            var currentEventProps = InitEvents(node, child.DistToParent);
            var currentPotential = LogPotential(child, fitnessMap, currentEventProps, ref thresholdAccept);
            // Now we perform the Metropolis-Hastings algorithm
            // and sample a set of events that give the closest agreement with
            // fitness given by SMITH
            for (int i = 0; i < _mcParams.NumSamplesTotal; i++)
            {
                // Reset the automatic acceptance
                thresholdAccept = false;
                var proposedEventProps = currentEventProps.ToList();
                // Select a random CNEventP to modify
                int index = _rnd.Next(proposedEventProps.Count);
                // Choose whether to swap the event entirely
                if (_rnd.NextDouble() < _mcParams.SwapEventP)
                {
                    // Get the new signature and the corresponding event
                    int mixtureIndex = Sampling.PickRandomIndex(_rnd, node.SigMixture);
                    var sig = _signatures[mixtureIndex];
                    SelectedSignatures[index] = sig;
                    var cnEventP = Sampling.PickRandomEventP(_rnd, sig.Events);
                    proposedEventProps[index] = Sampling.GenerateCNEventData(_rnd, node.Karyotype, cnEventP);
                }
                // Otherwise we modify some quantity of the event, but keep the event itself the same
                else
                {
                    // Keep the event type the same, but redo all parameters:
                    var cnEventP = proposedEventProps[index].EventP;
                    proposedEventProps[index] = Sampling.GenerateCNEventData(_rnd, node.Karyotype, cnEventP);

                }
                // With the newly selected event, we need to calculate the new
                // fitness of the clone
                var proposalPotential = LogPotential(child, fitnessMap, proposedEventProps, ref thresholdAccept);
                var logAcceptProb = proposalPotential.Item1 - currentPotential.Item1;
                if (logAcceptProb >= Math.Log(_rnd.NextDouble()))
                {
                    currentPotential = proposalPotential;
                    currentEventProps = proposedEventProps;
                    // Break out of the sampling if we have reached the threshold
                    // and have reached the minimum number of samples required
                    if (thresholdAccept && i > _mcParams.NumSamplesMin)
                        break;
                }
            }
            Console.WriteLine($"Final sampled fitness: {currentPotential.Item2}");
            //Console.WriteLine($"pot: {currentPotential.Item1}");

            // Finalize the mutated karyotype by applying the best-fit set of events
            // and move onto the next clone
            for (int mutNo = 0; mutNo < currentEventProps.Count; mutNo++)
            {
                Console.Write($"\rClone {counter}/{clones.Count - 1}. Event {mutNo + 1}/{child.DistToParent}.");
                var eventData = currentEventProps[mutNo];
                eventData.ApplyEvent(child.Karyotype);
                double newFitness = child.Karyotype.UpdateFitness(_geneLists, _fitness);
                int mutationCount = parentMutations + 1 + mutNo;
                double dFit = newFitness - oldFitness;
                var abberation = new CNEvent(child.CloneId, mutationCount, eventData.EventType, 
                    eventData.ToString(), dFit, newFitness);
                events.Add(abberation);
            }
            counter++;
            MCSampleCNEventsRec(child, clones, fitnessMap, events, ref counter);
        }
    }
}
