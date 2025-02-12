using SimChA.Computation;
using SimChA.Data;
using SimChA.EventData;
using SimChA.IO;

namespace SimChA.Simulation;

public class SASimulator : Simulator
{
    private EvoParams EvoParams { get; }

    public SASimulator(Random rnd, GenRef genRef, SimParams simParams, FitParams fitParams, EvoParams evoParams) 
        : base(rnd, genRef, simParams, fitParams)
    {
        EvoParams = evoParams;
    }
    
    private List<CNEventPars> GetEventPars(List<CNEventPars> pars, bool hasWGD)
    {
        if (EvoParams.EventCost <= 0 || !hasWGD)
        {
            return pars;
        }

        var newPars = new List<CNEventPars>(pars);
        foreach (var e in pars.Where(e => e.Type.ToString().EndsWith("Deletion")))
        {
            newPars[newPars.IndexOf(e)] = e with { Prob = e.Prob * EvoParams.EventCost };
        }
        var normalized = Factory.NormalizeEvents(newPars);
        return normalized;
    }

    private BaseEventData GetNewEvent(List<CNEventPars> cnEventPars, Karyotype kar, double currentFitness)
    {
        for (int tryNo = 0; tryNo <= EvoParams.MaxTries; tryNo++)
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
            if (Math.Exp(proposedFitness - currentFitness - EvoParams.Acceptance) >= Rnd.NextDouble())
            {
                return eventData;
            }
        }
        throw new Exception($"Could not generate a new event in {EvoParams.MaxTries} tries.");
    }

    bool MetBreakCondition(int time, int nEventsLeft)
        => !EvoParams.EvolveInTime ? nEventsLeft == 0 : time >= 1; 

    double GetMutRate(Karyotype kar)
        => EvoParams.EventRate * SampleStat.CalcPloidy(kar, GenRef);

    bool DidMutate(Karyotype kar)
        => !EvoParams.EvolveInTime || Rnd.NextDouble() < 1 - Math.Exp(-GetMutRate(kar));
    
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
        int distance = child.Distance > 0
            ? child.Distance
            : Sampling.SampleDistInt(Rnd, SimParams.RateDist, SimParams.RateMean);

        for (int mutNo = 1; mutNo <= distance; mutNo++)
        {
            Console.Write($"\rSample {child.CloneId}. Event {mutNo}/{distance}.".PadRight(80));
            var cnEventPars = GetEventPars(cnEventPs, hasWGD);
            var newEvent = GetNewEvent(cnEventPars, childKar, currentFit);
            var newKar = new Karyotype(childKar);
            newEvent.ApplyEvent(newKar);
            double proposedFit = newKar.UpdateFitness(GenRef, FitParams);
            double dFit = proposedFit - currentFit;
            var abberation = new CNEventDesc(newEvent.EventType, mutDepth + mutNo, newEvent.ToString(),
                dFit, proposedFit, mutNo);
            childEvs.Add(abberation);
                
            hasWGD |= newEvent.EventType == CNEventType.WholeGenomeDoubling;
            childKar = newKar;
            currentFit = proposedFit;
        }
        return (childKar, childEvs);
    }
}