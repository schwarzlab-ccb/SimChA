using SimChA.Computation;
using SimChA.Data;
using SimChA.EventData;
using SimChA.IO;

namespace SimChA.Simulation;

public class MatchSimulator(Random rnd, RefGen refGen, SimParams simParams, FitParams fitParams, EvoParams evoParams)
    : Simulator(rnd, refGen, simParams, fitParams)
{
    private EvoParams EvoParams { get; } = evoParams;

    private (Karyotype newKar, BaseEventData eventData, int numTries, string signature) GetNewEvent(
        List<CNEventPars> cnEventPars, Karyotype currentKar, double targetFitness, double decay)
    {
        double currentDist = Math.Abs(currentKar.FitnessVal - targetFitness);

        for (int tryNo = 0; tryNo <= EvoParams.MaxTries; tryNo++)
        {
            if (SampleStat.CalcPloidy(currentKar, RefGen) > 32)
                break;

            var cnEventP = Rnd.PickRndElem(cnEventPars);
            var eventData = Sampling.GenerateCNEventData(Rnd, currentKar, cnEventP);
            if (eventData == null)
                continue;

            var proposedKar = new Karyotype(currentKar);
            eventData.ApplyEvent(proposedKar);
            double proposedFitness = proposedKar.UpdateFitness(RefGen, FitParams);
            double proposedDist = Math.Abs(proposedFitness - targetFitness);

            // Accept if distance improves; decay makes acceptance stricter as fewer events remain
            double distImprovement = currentDist - proposedDist;
            if (Math.Exp(distImprovement - EvoParams.Acceptance - decay) > Rnd.NextDouble())
            {
                return (proposedKar, eventData, tryNo, cnEventP.Signature);
            }
        }
        return (currentKar, CreateSkipEvent(), EvoParams.MaxTries, "");
    }

    protected override (Karyotype childKar, List<CNEventDesc> childEvs) SampleEvents(
        Karyotype parentKar,
        CTreeNode cnChild,
        List<CNEventPars> cnEventPs,
        int mutDepth)
    {
        var currentKar = new Karyotype(parentKar);
        currentKar.UpdateFitness(RefGen, FitParams);
        var childEvs = new List<CNEventDesc>();
        double targetFit = cnChild.Fitness;

        int eventCount = SampleEventCount(cnChild);
        for (int evNo = 1; evNo <= eventCount; evNo++)
        {
            Console.Write($"\rSample {cnChild.CloneId}. Event {evNo}/{eventCount}.".PadRight(80));
            double oldFitness = currentKar.FitnessVal;
            // Decay increases as fewer events remain: 0 at start, Decay at last event
            double decay = eventCount > 1 
                ? EvoParams.Decay * (evNo - 1.0) / (eventCount - 1.0) 
                : 0;
            var (newKar, eventData, numTries, signature) = GetNewEvent(cnEventPs, currentKar, targetFit, decay);
            double newFit = newKar.FitnessVal;
            double dFit = newFit - oldFitness;
            var newEv = new CNEventDesc(eventData, mutDepth + evNo, dFit, newFit, numTries, Signature: signature);
            childEvs.Add(newEv);
            currentKar = newKar;
        }

        return (currentKar, childEvs);
    }
}