using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;
using EDists = Extreme.Statistics.Distributions;

namespace SimChA.Simulation;

public class MCSimulator : Simulator
{
    private FitnessParams Fitness { get; }
    private MCParams McParams { get; }
    private double TotalEventWeight {get; set;}
    public MCSimulator(
        Random rnd,
        GenRef genRef,
        FitnessParams fitnessParams, 
        MCParams mCParams) : base(rnd, genRef)
    {
        Fitness = fitnessParams;
        McParams = mCParams;
    }
    
    public override void SampleEvents(Sample sample)
    {
        if (sample.EventPars == null || !sample.EventPars.Any())
        {
            throw new Exception("No events to sample from.");
        }
        foreach( var e in sample.EventPars)
        {
            TotalEventWeight += e.Prob;
        }
        Counter = 1;
        var (root, childLoopUp) = CloneComp.CreateLookUp(sample.Clones);
        sample.Kars[root.CloneId] = new Karyotype(GenRef, sample.SexXX);
        ApplyCNEventsRec(sample, root, childLoopUp, 1);
    }
    
    public (double potential, bool accept) Potential(Karyotype kar, double targetFit, List<BaseEventData> events)
    {
        double eventPotentialTotal = 0.0;

        // Probability of picking each event and their corresponding signature
        foreach (var eventData in events)
        {
            eventData.ApplyEvent(kar);
            eventPotentialTotal += Math.Log(eventData.CNEventPars.Prob/TotalEventWeight);
        }
        double dFit = kar.UpdateFitness(GenRef, Fitness) - targetFit;
        // Variable to immediately quit the MC Sampling if we've reached enough accuracy
        bool accept = Math.Abs(dFit / targetFit) < McParams.ThresholdFit;

        double fitnessPotential = -McParams.ThetaFitness * Math.Abs(dFit);

        double potential = eventPotentialTotal + fitnessPotential;

        return (potential, accept);
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
    
    private List<BaseEventData> GetBestEvents(Sample sample, Karyotype kar, int nEvents, double targetFitness){
        // Generate a starting set of mutations and its potential
        var currentEvents = InitEvents(kar, nEvents, sample.EventPars);
        double currentPotential = Potential(new Karyotype(kar), targetFitness, currentEvents).potential;

        double bestPotential = currentPotential;
        var bestEvents = new List<BaseEventData>(currentEvents);
        for (int i = 0; i < McParams.NumSamplesTotal; i++)
        {
            var proposedEvents = GetNewProposal(sample, kar, currentEvents);
            // Calculate the new fitness of the proposed set of events on the clone
            (double proposalPotential, bool thresholdAccept) = Potential(new Karyotype(kar), targetFitness, proposedEvents);
            double acceptProb = proposalPotential - currentPotential;
            if (acceptProb >= Math.Log(Rnd.NextDouble()))
            {
                currentPotential = proposalPotential;
                currentEvents = proposedEvents;
                if (proposalPotential < bestPotential)
                {
                    bestPotential = proposalPotential;
                    bestEvents = proposedEvents;
                }
                // Break out of the sampling if we have reached the threshold
                // and have reached the minimum number of samples required
                if (thresholdAccept && i > McParams.NumSamplesMin)
                {
                    return currentEvents;
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
                
                var bestEvents = GetBestEvents(sample, childKar, child.Distance,child.FitnessTarget);

                for (int mutNo = 0; mutNo < bestEvents.Count; mutNo++)
                {
                    Console.Write($"\rSample {sample.SampleId}. Clone {Counter}/{clones.Count}. Event {mutNo + 1}/{child.Distance}.");
                    var eventData = bestEvents[mutNo];
                    eventData.ApplyEvent(childKar);
                    double newFitness = childKar.UpdateFitness(GenRef, Fitness);
                    double dFit = newFitness - oldFitness;
                    var abberation = new CNEventDesc(eventData.EventType, eventCount + mutNo, eventData.ToString(), dFit,
                        newFitness);
                    childEvs.Add(abberation);
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
