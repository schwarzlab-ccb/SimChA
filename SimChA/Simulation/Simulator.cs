using SimChA.DataTypes;
using SimChA.EventData;
using EDists = Extreme.Statistics.Distributions;

namespace SimChA.Simulation;

public class Simulator
{
    private readonly Random _rnd;
    private readonly FitnessParams _fitness;
    private readonly MCParams? _mcParams;
    private readonly Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> _geneLists;
    private int _counter;

    public Simulator(
        Random rnd,
        FitnessParams fitnessParams,
        MCParams mcParams,
        Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> geneLists)
    {
        _rnd = rnd;
        _fitness = fitnessParams;
        _mcParams = mcParams;
        _geneLists = geneLists;
    }
    
    public List<CNEvent> SampleEvents(Sample sample)
    {
        if (sample.EventPs == null || !sample.EventPs.Any())
        {
            throw new Exception("No events to sample from.");
        }
        var events = new List<CNEvent>();
        var (root, childLoopUp) = Clone.CreateLookUp(sample.Clones);
        _counter = 0;
        ApplyCNEventsRec(sample.EventPs, root, childLoopUp, events);
        return events;
    }

    private void ApplyCNEventsRec(List<CNEventP> eventsPs, Clone node, Dictionary<int, List<Clone>> clones, List<CNEvent> eventSeq)
    {
        foreach (var child in clones[node.CloneId])
        {
            child.Karyotype = node.CopyKaryotype();
            double oldFitness = node.Karyotype.FitnessVal;
            for (var mutNo = 0; mutNo < child.DistToParent; mutNo++)
            {
                Console.Write($"\rClone {_counter}/{clones.Count - 1}. Event {mutNo + 1}/{child.DistToParent}.");
                var eventP =Sampling.PickRndElem(_rnd, eventsPs);
                var eventData = Sampling.GenerateCNEventData(_rnd, child.Karyotype, eventP);
                eventData.ApplyEvent(child.Karyotype);
                double newFitness = child.Karyotype.UpdateFitness(_geneLists, _fitness);
                double dFit = newFitness - oldFitness;
                var abberation = new CNEvent(child.CloneId, eventP.Type, eventData.ToString(), dFit, newFitness);
                eventSeq.Add(abberation);
                oldFitness = newFitness;
            }
            _counter++;
            ApplyCNEventsRec(eventsPs, child, clones, eventSeq);
        }
    }
    
    public static List<Sample> SamplesFromProfiles(Dictionary<string, Karyotype> profiles)
    {
        var samples = new List<Sample>();
        foreach (var profile in profiles)
        {
            var clone = new Clone(0, -1, 0, profile.Value, 0);
            var sample = new Sample(profile.Key, profile.Value.SexXX, new List<Clone>{ clone });
            samples.Add(sample);
        }
        return samples;
    }

    public List<CNEvent> MCSampleEvents(Sample sample, Dictionary<int, double> fitnessMap)
    {        
        if (sample.EventPs == null || !sample.EventPs.Any())
        {
            throw new Exception("No events to sample from.");
        }
        List<CNEvent> events = new();
        var (root, childLoopUp) = Clone.CreateLookUp(sample.Clones);
        MCSampleCNEventsRec(sample.EventPs, root, childLoopUp, fitnessMap, events);
        return events;
    }

    public List<BaseEventData> InitEvents(Karyotype kar, int nMutations, List<CNEventP> cnEventPs)
    {
        var eventPs = Enumerable.Range(0, nMutations).Select(i => Sampling.PickRndElem(_rnd, cnEventPs));
        return eventPs.Select(e => Sampling.GenerateCNEventData(_rnd, kar, e)).ToList();
    }

    // The conditional (unnormalized) probability of this set of events occuring, 
    // given the individual events and the signature
    public double Potential(Clone node, Dictionary<int, double> fitnessMap,
        List<BaseEventData> events, ref bool thresholdAccept)
    {
        double eventPotentialTotal = 1.0;
        double targetFitness = fitnessMap[node.CloneId];

        // Create a dummy karyotype for the events to act on
        var karyotype = node.CopyKaryotype();

        // Probability of picking each event and their corresponding signature
        // (ignore normalization, divides out when we do the accept/reject)
        for (int i = 0; i < events.Count; i++)
        {
            events[i].ApplyEvent(karyotype);
            eventPotentialTotal *= events[i].EventP.Prob;
        }

        double newFitness = karyotype.UpdateFitness(_geneLists, _fitness);
        double dFit = newFitness - targetFitness;
        thresholdAccept = Math.Abs(dFit / targetFitness) < _mcParams.ThresholdFit;

        // Fitness potential is an exponential - exp[-theta * |fit - mean_fit|]
        double fitnessPotential = Math.Exp(-_mcParams.ThetaFitness * Math.Abs(dFit));
        // Gaussian form
        //double fitnessPotential = Math.Exp(-McParams.ThetaFitness*Math.Pow(dFit,2));

        return fitnessPotential * eventPotentialTotal;
    }
    
    private void MCSampleCNEventsRec(List<CNEventP> eventsPs, Clone node, Dictionary<int, List<Clone>> clones, 
        Dictionary<int, double> fitnessMap, List<CNEvent> events)
    {
        foreach (var child in clones[node.CloneId])
        {
            // Skip any clones that don't have any mutational distance from their
            // parent
            if (child.DistToParent == 0)
                continue;

            // Initialize all the relevant quantities
            double oldFitness = node.Karyotype.FitnessVal;

            // Parameters needed for the MH algorithm
            float alterEventStart = 0.5f;
            float alterEventLength = 0.5f;
            // The tracker to accept sets of events if they are within the threshold tolerance from
            // the target fitness
            bool thresholdAccept = false;

            // Generate a starting set of mutations and its potential
            var currentEventProps = InitEvents(node.Karyotype, child.DistToParent, eventsPs);
            double currentPotential = Potential(node, fitnessMap, currentEventProps, ref thresholdAccept);

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
                    var cnEventP = Sampling.PickRndElem(_rnd, eventsPs);
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
                double proposalPotential = Potential(node, fitnessMap, proposedEventProps, ref thresholdAccept);
                double acceptProb = proposalPotential / currentPotential;
                if (acceptProb >= _rnd.NextDouble())
                {
                    currentPotential = proposalPotential;
                    currentEventProps = proposedEventProps;
                    // Break out of the sampling if we have reached the threshold
                    // and have reached the minimum number of samples required
                    if (thresholdAccept && i > _mcParams.NumSamplesMin)
                        break;
                }
            }

            // Finalize the mutated karyotype by applying the best-fit set of events
            // and move onto the next clone
            // Copy its parent
            child.Karyotype = node.CopyKaryotype();
            for (int mutNo = 0; mutNo < currentEventProps.Count; mutNo++)
            {
                Console.Write($"\rClone {_counter}/{clones.Count - 1}. Event {mutNo + 1}/{child.DistToParent}.");
                var eventData = currentEventProps[mutNo];
                eventData.ApplyEvent(child.Karyotype);
                double newFitness = child.Karyotype.UpdateFitness(_geneLists, _fitness);
                double dFit = newFitness - oldFitness;
                var abberation = new CNEvent(child.CloneId, eventData.EventType, eventData.ToString(), dFit, newFitness);
                events.Add(abberation);
            }
            _counter++;
            MCSampleCNEventsRec(eventsPs, child, clones, fitnessMap, events);
        }
    }
}