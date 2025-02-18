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

    private (Karyotype newKar, BaseEventData eventData) GetNewEvent(List<CNEventPars> cnEventPars, Karyotype currentKar)
    {
        for (int tryNo = 0; tryNo <= SAParams.MaxTries; tryNo++)
        {
            var cnEventP = Rnd.PickRndElem(cnEventPars);
            var eventData = Sampling.GenerateCNEventData(Rnd, currentKar, cnEventP);
            if (eventData != null)
            {
                var proposedKar = new Karyotype(currentKar);
                eventData.ApplyEvent(proposedKar);
                double proposedFitness = proposedKar.UpdateFitness(GenRef, FitParams);
                if (Math.Exp(proposedFitness - currentKar.FitnessVal - SAParams.Acceptance) > Rnd.NextDouble())
                {
                    return (proposedKar, eventData);
                }
            }
        }
        return (currentKar, CreatePassEvent());
    }
    
    protected override (Karyotype childKar, List<CNEventDesc> childEvs) SampleEvents(
        Karyotype currentKar, 
        CTreeNode cnChild, 
        List<CNEventPars> cnEventPs, 
        int mutDepth)
    {
        var childEvs = new List<CNEventDesc>();
        bool hasWGD = false;
        int eventCount = SampleEventCount(cnChild);
        for (int evNo = 1; evNo <= eventCount; evNo++)
        {
            Console.Write($"\rSample {cnChild.CloneId}. Event {evNo}/{eventCount}.".PadRight(80));
            var cnEventPars = GetEventPars(cnEventPs, hasWGD);
            var (newKar, eventData) = GetNewEvent(cnEventPars, currentKar);
            double newFit = newKar.UpdateFitness(GenRef, FitParams);
            double dFit = newFit - currentKar.FitnessVal;
            var newEv = new CNEventDesc(eventData, mutDepth + evNo, dFit, newFit);
            
            childEvs.Add(newEv);
            hasWGD |= eventData.EventType == CNEventType.WholeGenomeDoubling;
            currentKar = newKar;
        }
        return (currentKar, childEvs);
    }
}