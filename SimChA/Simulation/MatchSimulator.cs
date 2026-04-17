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
        List<CNEventPars> cnEventPars, Karyotype currentKar, double targetFitness)
    {
        Karyotype? bestKar = null;
        BaseEventData? bestEvent = null;
        string bestSignature = "";
        double bestDist = double.PositiveInfinity;
        int bestTry = 0;

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

            if (proposedDist < bestDist)
            {
                bestKar = proposedKar;
                bestEvent = eventData;
                bestSignature = cnEventP.Signature;
                bestDist = proposedDist;
                bestTry = tryNo;
            }
        }
        return bestKar is not null 
            ? (bestKar, bestEvent!, bestTry, bestSignature) 
            : (currentKar, CreateSkipEvent(), EvoParams.MaxTries, "");
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
            var (newKar, eventData, numTries, signature) = GetNewEvent(cnEventPs, currentKar, targetFit);
            double newFit = newKar.FitnessVal;
            double dFit = newFit - oldFitness;
            var newEv = new CNEventDesc(eventData, mutDepth + evNo, dFit, newFit, numTries, Signature: signature);
            childEvs.Add(newEv);
            currentKar = newKar;
        }

        return (currentKar, childEvs);
    }
}