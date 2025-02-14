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
    {
        SAParams = saParams;
    }
    
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
        throw new Exception($"Could not generate a new event in {SAParams.MaxTries} tries.");
    }

    bool MetBreakCondition(int time, int nEventsLeft)
        => !SAParams.EvolveInTime ? nEventsLeft == 0 : time >= 1; 

    double GetMutRate(Karyotype kar)
        => SAParams.EventRate * SampleStat.CalcPloidy(kar, GenRef);

    bool DidMutate(Karyotype kar)
        => !SAParams.EvolveInTime || Rnd.NextDouble() < 1 - Math.Exp(-GetMutRate(kar));
    
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

        double rate = SampleRate(child);
        double time = 0.0;
        do
        {
            time += rate;
            int evNo = childEvs.Count + 1;
            Console.Write($"\rSample {child.CloneId}. Event {evNo}/{child.Distance}.".PadRight(80));
            var cnEventPars = GetEventPars(cnEventPs, hasWGD);
            var newEvent = GetNewEvent(cnEventPars, childKar, currentFit);
            var newKar = new Karyotype(childKar);
            newEvent.ApplyEvent(newKar);
            double proposedFit = newKar.UpdateFitness(GenRef, FitParams);
            double dFit = proposedFit - currentFit;
            var abberation = new CNEventDesc(newEvent.EventType, mutDepth + evNo, newEvent.ToString(),
                dFit, proposedFit, time);
            
            childEvs.Add(abberation);
                
            hasWGD |= newEvent.EventType == CNEventType.WholeGenomeDoubling;
            childKar = newKar;
            currentFit = proposedFit;
        } while (time <= 1 - double.Epsilon);
        return (childKar, childEvs);
    }
}