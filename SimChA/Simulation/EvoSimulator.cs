using SimChA.Computation;
using SimChA.Data;
using SimChA.EventData;
using SimChA.IO;

namespace SimChA.Simulation;

public class EvoSimulator(Random rnd, GenRef genRef, SimParams simParams, FitParams fitParams, EvoParams evoParams)
    : Simulator(rnd, genRef, simParams, fitParams)
{
    private EvoParams EvoParams { get; } = evoParams;

    private (Karyotype newKar, BaseEventData eventData, int numTries) GetNewEvent(List<CNEventPars> cnEventPars, Karyotype currentKar)
    {
        for (int tryNo = 0; tryNo <= EvoParams.MaxTries; tryNo++)
        {
            if (SampleStat.CalcPloidy(currentKar, GenRef) > 32)
            {
                break;
            }
            var cnEventP = Rnd.PickRndElem(cnEventPars);
            var eventData = Sampling.GenerateCNEventData(Rnd, currentKar, cnEventP);
            if (eventData != null)
            {
                var proposedKar = new Karyotype(currentKar);
                eventData.ApplyEvent(proposedKar);
                double proposedFitness = proposedKar.UpdateFitness(GenRef, FitParams);
                if (Math.Exp(proposedFitness - currentKar.FitnessVal - EvoParams.Acceptance) > Rnd.NextDouble())
                {
                    return (proposedKar, eventData, tryNo);
                }
            }
        }
        return (currentKar, CreatePassEvent(), EvoParams.MaxTries);
    }
    
    protected override (Karyotype childKar, List<CNEventDesc> childEvs) SampleEvents(
        Karyotype currentKar, 
        CTreeNode cnChild, 
        List<CNEventPars> cnEventPs, 
        int mutDepth)
    {
        var childEvs = new List<CNEventDesc>();
        int eventCount = SampleEventCount(cnChild);
        for (int evNo = 1; evNo <= eventCount; evNo++)
        {
            Console.Write($"\rSample {cnChild.CloneId}. Event {evNo}/{eventCount}.".PadRight(80));
            var (newKar, eventData, numTries) = GetNewEvent(cnEventPs, currentKar);
            double newFit = newKar.FitnessVal;
            double dFit = newFit - currentKar.FitnessVal;
            var newEv = new CNEventDesc(eventData, mutDepth + evNo, dFit, newFit, numTries);
            
            childEvs.Add(newEv);
            currentKar = newKar;
        }
        return (currentKar, childEvs);
    }
}