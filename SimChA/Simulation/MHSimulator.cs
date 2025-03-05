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

    private double Gaussian(double f, double ft)
        => Math.Exp(-Math.Pow(f - ft, 2.0)/(2*MHParams.ThetaFitness));

    private double AcceptProb(double newFit, double oldFit, double targetFitness)
        => Gaussian(newFit, targetFitness) / Gaussian(oldFit, targetFitness);

    private double AcceptProb(double newFit, double oldFit)
        => Math.Min(1, Math.Exp((newFit - oldFit)/MHParams.ThetaFitness));

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
        int index = Rnd.Next(proposedEvents.Count);
        var newData = Sampling.GenerateCNEventData(Rnd, kar, Rnd.PickRndElem(cnEventPs));
        proposedEvents[index] = newData ?? throw new Exception("Failed to generate new event data.");
        return proposedEvents;
    }

    private List<BaseEventData> GenEventsForTargetFitness(List<CNEventPars> cnEventPs, Karyotype kar, int nEvents,
        double targetFitness)
    {
        // Generate a starting set of mutations and its potential
        var currentEvents = InitEvents(kar, nEvents, cnEventPs);
        double currentFitness = GetFitness(new Karyotype(kar), currentEvents);
        double bestDiff = double.PositiveInfinity;
        var bestEvents = new List<BaseEventData>(currentEvents);

        for (int i = 0; i < MHParams.NumIterations; i++)
        {
            var proposedEvents = GetNewProposal(cnEventPs, kar, currentEvents);
            double proposedFitness = GetFitness(new Karyotype(kar), proposedEvents);
            if (AcceptProb(proposedFitness, currentFitness, targetFitness) > Rnd.NextDouble())
            {
                currentEvents = proposedEvents;
                double proposedDiff = Math.Abs(proposedFitness - targetFitness);
                if (proposedDiff < bestDiff)
                {
                    bestDiff = proposedDiff;
                    bestEvents = proposedEvents;
                }
            }
        }

        return bestEvents;
    }

    private List<BaseEventData> GenEventsForMaxFitness(List<CNEventPars> cnEventPs, Karyotype kar, int nEvents)
    {
        // Generate a starting set of mutations and its potential
        var currentEvents = InitEvents(kar, nEvents, cnEventPs);
        double currentFitness = GetFitness(new Karyotype(kar), currentEvents);
        double bestFitness = currentFitness;
        var bestEvents = new List<BaseEventData>(currentEvents);

        for (int i = 0; i < MHParams.NumIterations; i++)
        {
            var proposedEvents = GetNewProposal(cnEventPs, kar, currentEvents);
            var proposedFitness = GetFitness(new Karyotype(kar), proposedEvents);
            if (AcceptProb(proposedFitness, currentFitness) > Rnd.NextDouble())
            {
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
        CTreeNode cnChild,
        List<CNEventPars> cnEventPs,
        int mutDepth)
    {
        var childKar = new Karyotype(parentKar);
        var childEvs = new List<CNEventDesc>();
        double targetFit = SampleFit(cnChild);
        double oldFitness = childKar.FitnessVal;

        int eventCount = SampleEventCount(cnChild);
        var bestEvents = MHParams.MatchFitness
            ? GenEventsForTargetFitness(cnEventPs, childKar, eventCount, targetFit)
            : GenEventsForMaxFitness(cnEventPs, childKar, eventCount);

        for (int evNo = 1; evNo <= eventCount; evNo++)
        {
            Console.Write($"\rSample {cnChild.CloneId}. Event {evNo}/{eventCount}.".PadRight(80));

            var eventData = bestEvents[evNo - 1];
            eventData.ApplyEvent(childKar);
            double newFitness = childKar.UpdateFitness(GenRef, FitParams);
            double dFit = newFitness - oldFitness;
            var newEv = new CNEventDesc(eventData, mutDepth + evNo, dFit, newFitness);
            childEvs.Add(newEv);
            oldFitness = newFitness;
        }

        return (childKar, childEvs);
    }
}