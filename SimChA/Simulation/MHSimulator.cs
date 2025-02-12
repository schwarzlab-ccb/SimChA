using SimChA.Computation;
using SimChA.Data;
using SimChA.EventData;
using SimChA.IO;

namespace SimChA.Simulation;

public class MHSimulator : Simulator
{
    private MHParams MHParams { get; }

    public MHSimulator(Random rnd, GenRef genRef, SimParams simParams, FitParams fitParams, MHParams mhParams)
        : base(rnd, genRef, simParams, fitParams)
    {
        MHParams = mhParams;
    }

    private List<BaseEventData> InitEvents(Karyotype kar, int nMutations, List<CNEventPars> cnEventPs)
    {
        var eventPs = Enumerable.Range(0, nMutations).Select(_ => Rnd.PickRndElem(cnEventPs));
        return eventPs.Select(
            e => Sampling.GenerateCNEventData(Rnd, kar, e)
                 ?? throw new Exception($"Failed to generate event data for {e}.")
        ).ToList();
    }

    private double GetFitnessPotential(double fitness, double targetFitness)
        => -MHParams.ThetaFitness * Math.Abs((fitness - targetFitness) / targetFitness);

    private double CalculatePotential(double proposedFitness)
        => MHParams.ThetaFitness * proposedFitness;

    private double CalculatePotential(double proposedFitness, double targetFitness)
        => GetFitnessPotential(proposedFitness, targetFitness);

    private double GetFitness(Karyotype kar, List<BaseEventData> events)
    {
        // Probability of picking each event and their corresponding signature
        foreach (var eventData in events)
        {
            eventData.ApplyEvent(kar);
        }

        return kar.UpdateFitness(GenRef, FitParams);
    }

    private List<BaseEventData> GetNewProposal(List<CNEventPars> cnEventPs, Karyotype kar,
        List<BaseEventData> oldEvents)
    {
        var proposedEvents = oldEvents.ToList();
        // Select a random CNEventPars to modify
        int index = Rnd.Next(proposedEvents.Count);
        var cnEventP = Rnd.NextDouble() < MHParams.SwapEventP
            ? Rnd.PickRndElem(cnEventPs)
            : proposedEvents[index].CNEventPars;
        var newData = Sampling.GenerateCNEventData(Rnd, kar, cnEventP);
        proposedEvents[index] = newData ?? throw new Exception("Failed to generate new event data.");
        return proposedEvents;
    }

    private List<BaseEventData> GenEventsForTargetFitness(List<CNEventPars> cnEventPs, Karyotype kar, int nEvents,
        double targetFitness)
    {
        // Generate a starting set of mutations and its potential
        var currentEvents = InitEvents(kar, nEvents, cnEventPs);
        double currentFitness = GetFitness(new Karyotype(kar), currentEvents);
        double currentPotential = CalculatePotential(currentFitness, targetFitness);
        double bestDiff = 1000.0;
        var bestEvents = new List<BaseEventData>(currentEvents);

        for (int i = 0; i < MHParams.NumSamplesTotal; i++)
        {
            var proposedEvents = GetNewProposal(cnEventPs, kar, currentEvents);
            double fitness = GetFitness(new Karyotype(kar), proposedEvents);
            bool thresholdAccept = Math.Abs(1.0 - fitness / targetFitness) < MHParams.ThresholdFit;
            // Calculate the new fitness of the proposed set of events on the clone
            double proposalPotential = CalculatePotential(currentFitness, targetFitness);
            double acceptProb = proposalPotential - currentPotential;
            if (acceptProb >= Math.Log(Rnd.NextDouble()))
            {
                currentPotential = proposalPotential;
                currentEvents = proposedEvents;
                double proposedDiff = Math.Abs(fitness - targetFitness);
                if (proposedDiff < bestDiff)
                {
                    bestDiff = proposedDiff;
                    bestEvents = proposedEvents;
                }

                // Break out of the sampling if we have reached the threshold
                // and have reached the minimum number of samples required
                if (thresholdAccept && i > MHParams.NumSamplesMin) break;
            }
        }

        return bestEvents;
    }

    private List<BaseEventData> GenEventsForMaxFitness(List<CNEventPars> cnEventPs, Karyotype kar, int nEvents)
    {
        // Generate a starting set of mutations and its potential
        var currentEvents = InitEvents(kar, nEvents, cnEventPs);
        double currentFitness = GetFitness(new Karyotype(kar), currentEvents);
        double currentPotential = CalculatePotential(currentFitness);
        double bestFitness = currentFitness;
        var bestEvents = new List<BaseEventData>(currentEvents);

        for (int i = 0; i < MHParams.NumSamplesTotal; i++)
        {
            var proposedEvents = GetNewProposal(cnEventPs, kar, currentEvents);
            var proposedFitness = GetFitness(new Karyotype(kar), proposedEvents);
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

    protected override (Karyotype childKar, List<CNEventDesc> childEvs) SampleEvents(
        Karyotype parentKar,
        CTreeNode child,
        List<CNEventPars> cnEventPs,
        int mutDepth)
    {
        var childKar = new Karyotype(parentKar);
        var childEvs = new List<CNEventDesc>();
        int distance = SampleDist(child);
        double targetFit = SampleFit(child);
        double oldFitness = childKar.FitnessVal;

        var bestEvents = MHParams.MatchFitness
            ? GenEventsForTargetFitness(cnEventPs, childKar, distance, targetFit)
            : GenEventsForMaxFitness(cnEventPs, childKar, distance);

        for (int mutNo = 0; mutNo < bestEvents.Count; mutNo++)
        {
            Console.Write($"\rSample {child.CloneId}. Event {mutNo}/{child.Distance}.".PadRight(80));

            var eventData = bestEvents[mutNo];
            eventData.ApplyEvent(childKar);
            double newFitness = childKar.UpdateFitness(GenRef, FitParams);
            double dFit = newFitness - oldFitness;
            var newEv = new CNEventDesc(eventData.EventType, mutDepth + mutNo, eventData.ToString(),
                dFit, newFitness, mutNo);
            childEvs.Add(newEv);
            oldFitness = newFitness;
        }

        return (childKar, childEvs);
    }
}