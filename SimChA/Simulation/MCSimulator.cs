using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;
using EDists = Extreme.Statistics.Distributions;

namespace SimChA.Simulation;

public class MCSimulator : Simulator
{
    private FitnessParams FitnessParams { get; }
    private MCParams McParams { get; }
    public MCSimulator(
        Random rnd,
        GenRef genRef,
        FitnessParams fitnessParams, 
        MCParams mCParams) : base(rnd, genRef)
    {
        FitnessParams = fitnessParams;
        McParams = mCParams;
    }
    
    public override void SampleEvents(Sample sample)
    {
        if (sample.EventPars == null || !sample.EventPars.Any())
        {
            throw new Exception("No events to sample from.");
        }
        Counter = 1;
        var (root, childLoopUp) = CloneComp.CreateLookUp(sample.Clones);
        sample.Kars[root.CloneId] = new Karyotype(GenRef, sample.Sex);
        ApplyCNEventsRec(sample, root, childLoopUp, 1);
    }
    
    public (double potential, bool accept, double fitness) Potential(Karyotype kar, double targetFit, List<BaseEventData> events)
    {
        double eventPotentialTotal = 0.0;

        // Probability of picking each event and their corresponding signature
        foreach (var eventData in events)
        {
            eventData.ApplyEvent(kar);
            if (eventData.CNEventPars.Size > 0 && McParams.IncludeSize)
            {
                eventPotentialTotal += Math.Log(eventData.GetProb());
            }
            eventPotentialTotal += Math.Log(eventData.CNEventPars.Prob);
        }
        double newFitness = kar.UpdateFitness(GenRef, FitnessParams);
        double dFit = newFitness - targetFit;
        // Variable to immediately quit the MC Sampling if we've reached enough accuracy
        bool accept = Math.Abs(dFit / targetFit) < McParams.ThresholdFit;

        double fitnessPotential = -McParams.ThetaFitness * Math.Abs(dFit)/targetFit;

        double potential = fitnessPotential;
        if (McParams.IncludeProb) {
            potential += eventPotentialTotal;
        }

        return (potential, accept, newFitness);
    }

    public double CalculatePotential(double proposedFitness, double targetFitness, List<BaseEventData> events)
        => McParams.IncludeProb ? GetEventPotential(events) + GetFitnessPotential(proposedFitness, targetFitness) : GetFitnessPotential(proposedFitness, targetFitness);
    
    public double GetFitnessPotential(double fitness, double targetFitness)
    {
        double dFit = fitness - targetFitness;
        return -McParams.ThetaFitness * Math.Abs(dFit/targetFitness);
    }

    public double GetEventPotential(List<BaseEventData> events)
    {
        double eventPotentialTotal = 0.0;
        // Probability of picking each event and their corresponding signature
        foreach (var eventData in events)
        {
            if (eventData.CNEventPars.Size > 0 && McParams.IncludeSize)
            {
                eventPotentialTotal += Math.Log(eventData.GetProb());
            }
            eventPotentialTotal += Math.Log(eventData.CNEventPars.Prob);
        }
        return eventPotentialTotal;
    }

    public double CalculatePotential(double proposedFitness, List<BaseEventData> events)
    {
        double fitnessPotential = -McParams.ThetaFitness * proposedFitness;

        return McParams.IncludeProb ? GetEventPotential(events) + fitnessPotential : fitnessPotential;
    }

    public double GetFitness(Karyotype kar, List<BaseEventData> events)
    {
        // Probability of picking each event and their corresponding signature
        foreach (var eventData in events)
        {
            eventData.ApplyEvent(kar);
        }
        return kar.UpdateFitness(GenRef, FitnessParams);
    }

    private List<BaseEventData> GetNewProposal(Sample sample, Karyotype kar, List<BaseEventData> oldEvents)
    {
        var proposedEvents = oldEvents.ToList();
        // Select a random CNEventPars to modify
        int index = Rnd.Next(proposedEvents.Count);
        var cnEventP = Rnd.NextDouble() < McParams.SwapEventP 
            ? Rnd.PickRndElem(sample.EventPars) 
            : proposedEvents[index].CNEventPars;
        var newData = Sampling.GenerateCNEventData(Rnd, kar, cnEventP);
        if (newData != null)
        {
            proposedEvents[index] = newData;
        }
        else
        {
            Console.WriteLine("Failed to generate new event data.");
        }
        return proposedEvents;
    }
    
    private List<BaseEventData> GetEventsFromTargetFitness(Sample sample, Karyotype kar, int nEvents, double targetFitness)
    {
        // Generate a starting set of mutations and its potential
        var currentEvents = InitEvents(kar, nEvents, sample.EventPars);
        var currentFitness = GetFitness(new Karyotype(kar), currentEvents);
        double currentPotential = CalculatePotential(currentFitness, targetFitness, currentEvents);
        var best_diff = 1000.0;
        var bestEvents = new List<BaseEventData>(currentEvents);

        for (int i = 0; i < McParams.NumSamplesTotal; i++)
        {
            var proposedEvents = GetNewProposal(sample, kar, currentEvents);
            var fitness = GetFitness(new Karyotype(kar), proposedEvents);
            var thresholdAccept = Math.Abs(1.0 - fitness/targetFitness) < McParams.ThresholdFit;
            // Calculate the new fitness of the proposed set of events on the clone
            double proposalPotential = CalculatePotential(currentFitness, targetFitness, proposedEvents);
            double acceptProb = proposalPotential - currentPotential;
            if (acceptProb >= Math.Log(Rnd.NextDouble()))
            {
                currentPotential = proposalPotential;
                currentEvents = proposedEvents;
                var proposed_diff = Math.Abs(fitness - targetFitness);
                if (proposed_diff < best_diff)
                {
                    best_diff = proposed_diff;
                    bestEvents = proposedEvents;
                }
                // Break out of the sampling if we have reached the threshold
                // and have reached the minimum number of samples required
                if (thresholdAccept && i > McParams.NumSamplesMin) break;
            }
        }
        return bestEvents;
    }

    private List<BaseEventData> GetEventsFromMaximumFitness(Sample sample, Karyotype kar, int nEvents)
    {
        // Generate a starting set of mutations and its potential
        var currentEvents = InitEvents(kar, nEvents, sample.EventPars);
        var currentFitness = GetFitness(new Karyotype(kar), currentEvents);
        double currentPotential = CalculatePotential(currentFitness, currentEvents);
        var bestFitness = currentFitness;
        var bestEvents = new List<BaseEventData>(currentEvents);

        for (int i = 0; i < McParams.NumSamplesTotal; i++)
        {
            var proposedEvents = GetNewProposal(sample, kar, currentEvents);
            var fitness = GetFitness(new Karyotype(kar), proposedEvents);
            // Calculate the new fitness of the proposed set of events on the clone
            double proposalPotential = CalculatePotential(currentFitness, proposedEvents);
            double acceptProb = proposalPotential - currentPotential;
            if (acceptProb >= Math.Log(Rnd.NextDouble()))
            {
                currentPotential = proposalPotential;
                currentEvents = proposedEvents;
                if (currentFitness > bestFitness)
                {
                    bestFitness = currentFitness;
                    bestEvents = proposedEvents;
                }
            }
        }
        return bestEvents;
    }

    private void ApplyCNEventsRec(Sample sample, CloneIn node, IReadOnlyDictionary<int, 
        List<CloneIn>> clones, int eventCount)
    {
        foreach (var child in clones[node.CloneId])
        {
            var childKar = new Karyotype(sample.Kars[node.CloneId]);
            sample.Kars[child.CloneId] = childKar;
            var childEvs = new List<CNEventDesc>();
            sample.EventDescs[child.CloneId] = childEvs;
            
            if (child.Distance > 0)
            {
                double oldFitness = childKar.FitnessVal;
                
                var bestEvents = GetEventsFromTargetFitness(sample, childKar, child.Distance, child.FitnessTarget);

                for (int mutNo = 0; mutNo < bestEvents.Count; mutNo++)
                {
                    Console.Write($"\rSample {sample.SampleId}. Clone {Counter}/{clones.Count}. Event {mutNo + 1}/{child.Distance}.");
                    var eventData = bestEvents[mutNo];
                    eventData.ApplyEvent(childKar);
                    double newFitness = childKar.UpdateFitness(GenRef, FitnessParams);
                    double dFit = newFitness - oldFitness;
                    var abberation = new CNEventDesc(eventData.EventType, eventCount + mutNo, eventData.ToString(), dFit,
                        newFitness);
                    childEvs.Add(abberation);
                    oldFitness = newFitness;
                }
                Counter++;
                if (child.CloneId != node.CloneId)
                {
                    ApplyCNEventsRec(sample, child, clones, eventCount + child.Distance);
                }
            }
        }
    }
}
