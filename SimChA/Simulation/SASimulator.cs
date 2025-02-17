using SimChA.Computation;
using SimChA.Data;
using SimChA.EventData;
using SimChA.IO;

namespace SimChA.Simulation;

public class SASimulator : Simulator
{
    private SAParams SAParams { get; }

    public SASimulator(Random rnd, GenRef genRef, SimParams simParams, FitParams fitParams, SAParams saParams) 
        : base(rnd, genRef, simParams, fitParams) 
        => SAParams = saParams;

    private List<CNEventPars> GetEventPars(List<CNEventPars> pars, bool hasWGD)
    {
        if (SAParams.EventCost <= 0 || !hasWGD)
        {
            return pars;
        }

        var newPars = new List<CNEventPars>(pars);
        foreach (var e in pars.Where(e => e.Type.ToString().EndsWith("Deletion")))
        {
            newPars[newPars.IndexOf(e)] = e with { Prob = e.Prob * SAParams.EventCost };
        }
        var normalized = Factory.NormalizeEvents(newPars);
        return normalized;
    }

    private BaseEventData GetNewEvent(List<CNEventPars> cnEventPars, Karyotype kar, double currentFitness)
    {
        for (int tryNo = 0; tryNo <= SAParams.MaxTries; tryNo++)
        {
            var cnEventP = Rnd.PickRndElem(cnEventPars);
            var eventData = Sampling.GenerateCNEventData(Rnd, kar, cnEventP);
            if (eventData == null)
            {
                continue;
            }
            var proposedKar = new Karyotype(kar);
            eventData.ApplyEvent(proposedKar);
            double proposedFitness = proposedKar.UpdateFitness(GenRef, FitParams);
            if (Math.Exp(proposedFitness - currentFitness - SAParams.Acceptance) >= Rnd.NextDouble())
            {
                return eventData;
            }
        }
        return CreatePassEvent();
    }
    
    protected override (Karyotype childKar, List<CNEventDesc> childEvs) SampleEvents(
        Karyotype parentKar, 
        CTreeNode child, 
        List<CNEventPars> cnEventPs, 
        int mutDepth)
    {
        var childKar = new Karyotype(parentKar);
        var childEvs = new List<CNEventDesc>();
        double currentFit = childKar.FitnessVal;
        bool hasWGD = false;
        int eventCount = SampleEventCount(child);
        for (int evNo = 1; evNo <= eventCount; evNo++)
        {
            Console.Write($"\rSample {child.CloneId}. Event {evNo}/{eventCount}".PadRight(80));
            var cnEventPars = GetEventPars(cnEventPs, hasWGD);
            var eventData = GetNewEvent(cnEventPars, childKar, currentFit);
            var newKar = new Karyotype(childKar);
            eventData.ApplyEvent(newKar);
            double proposedFit = newKar.UpdateFitness(GenRef, FitParams);
            double dFit = proposedFit - currentFit;
            var newEv = new CNEventDesc(eventData, mutDepth + evNo, dFit, proposedFit);
            
            childEvs.Add(newEv);
            hasWGD |= eventData.EventType == CNEventType.WholeGenomeDoubling;
            childKar = newKar;
            currentFit = proposedFit;
        }
        return (childKar, childEvs);
    }
}