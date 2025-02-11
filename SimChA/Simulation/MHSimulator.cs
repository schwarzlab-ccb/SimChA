using SimChA.Computation;
using SimChA.Data;
using SimChA.EventData;
using SimChA.IO;

namespace SimChA.Simulation;

public class MHSimulator : Simulator
{
    private MHParams MHParams { get; }
    
    public MHSimulator(
        Random rnd,
        GenRef genRef,
        SimParams simParams,
        FitParams fitParams, 
        MHParams mhParams) : base(rnd, genRef, simParams, fitParams)
    {
        MHParams = mhParams;
    }

    public List<BaseEventData> InitEvents(Karyotype kar, int nMutations, List<CNEventPars> cnEventPs)
    {
        var eventPs = Enumerable.Range(0, nMutations).Select(_ => Rnd.PickRndElem(cnEventPs));
        return eventPs.Select(
            e => Sampling.GenerateCNEventData(Rnd, kar, e) ?? throw new Exception($"Failed to generate event data for {e}.")
        ).ToList();
    }
    
    public double GetFitnessPotential(double fitness, double targetFitness)
    {
        double dFit = fitness - targetFitness;
        return -MHParams.ThetaFitness * Math.Abs(dFit/targetFitness);
    }
    
    public double CalculatePotential(double proposedFitness)
        => MHParams.ThetaFitness * proposedFitness;
    
    public double CalculatePotential(double proposedFitness, double targetFitness)
        => GetFitnessPotential(proposedFitness, targetFitness);

    public double GetFitness(Karyotype kar, List<BaseEventData> events)
    {
        // Probability of picking each event and their corresponding signature
        foreach (var eventData in events)
            eventData.ApplyEvent(kar);
        return kar.UpdateFitness(GenRef, FitParams);
    }

    private List<BaseEventData> GetNewProposal(Sample sample, Karyotype kar, List<BaseEventData> oldEvents)
    {
        var proposedEvents = oldEvents.ToList();
        // Select a random CNEventPars to modify
        int index = Rnd.Next(proposedEvents.Count);
        var cnEventP = Rnd.NextDouble() < MHParams.SwapEventP 
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
    
    private List<BaseEventData> GenEventsForTargetFitness(Sample sample, Karyotype kar, int nEvents, double targetFitness)
    {
        // Generate a starting set of mutations and its potential
        var currentEvents = InitEvents(kar, nEvents, sample.EventPars);
        var currentFitness = GetFitness(new Karyotype(kar), currentEvents);
        double currentPotential = CalculatePotential(currentFitness, targetFitness);
        var best_diff = 1000.0;
        var bestEvents = new List<BaseEventData>(currentEvents);

        for (int i = 0; i < MHParams.NumSamplesTotal; i++)
        {
            var proposedEvents = GetNewProposal(sample, kar, currentEvents);
            var fitness = GetFitness(new Karyotype(kar), proposedEvents);
            var thresholdAccept = Math.Abs(1.0 - fitness/targetFitness) < MHParams.ThresholdFit;
            // Calculate the new fitness of the proposed set of events on the clone
            double proposalPotential = CalculatePotential(currentFitness, targetFitness);
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
                if (thresholdAccept && i > MHParams.NumSamplesMin) break;
            }
        }
        return bestEvents;
    }

    private List<BaseEventData> GenEventsForMaxFitness(Sample sample, Karyotype kar, int nEvents)
    {
        // Generate a starting set of mutations and its potential
        var currentEvents = InitEvents(kar, nEvents, sample.EventPars);
        var currentFitness = GetFitness(new Karyotype(kar), currentEvents);
        double currentPotential = CalculatePotential(currentFitness);
        var bestFitness = currentFitness;
        var bestEvents = new List<BaseEventData>(currentEvents);

        var fitList = new List<double>{currentFitness};

        for (int i = 0; i < MHParams.NumSamplesTotal; i++)
        {
            var proposedEvents = GetNewProposal(sample, kar, currentEvents);
            var proposedFitness = GetFitness(new Karyotype(kar), proposedEvents);
            fitList.Add(proposedFitness);
            // Calculate the new fitness of the proposed set of events on the clone
            double proposalPotential = CalculatePotential(proposedFitness);
            double acceptProb = proposalPotential - currentPotential;
            if (acceptProb >= Math.Log(Rnd.NextDouble()))
            {
                currentPotential = proposalPotential;
                currentEvents = proposedEvents;
                if (proposedFitness > bestFitness)
                {
                    bestFitness = proposedFitness;
                    bestEvents = proposedEvents;
                }
            }
        }
        return bestEvents;
    }

    protected virtual void ApplyCNEventsRec(
        CTreeNode parent, 
        List<CTreeNode> cloneTree, 
        List<CNEventPars> cnEventPs,
        Dictionary<string, double> mixture,
        List<Sample> sampleList,
        Karyotype parentKar,
        int mutDepth)
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
                
                var bestEvents = MHParams.MatchFitness
                    ? GenEventsForTargetFitness(sample, childKar, child.Distance, child.Fitness)
                    : GenEventsForMaxFitness(sample, childKar, child.Distance);

                for (int mutNo = 0; mutNo < bestEvents.Count; mutNo++)
                {
                    Console.Write($"\rSample {sample.SampleId}. Clone {Counter}/{clones.Count}. Event {mutNo + 1}/{child.Distance}.");
                    var eventData = bestEvents[mutNo];
                    eventData.ApplyEvent(childKar);
                    double newFitness = childKar.UpdateFitness(GenRef, FitParams);
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
