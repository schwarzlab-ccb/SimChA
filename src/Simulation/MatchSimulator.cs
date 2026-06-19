using SimChA.Computation;
using SimChA.Data;
using SimChA.EventData;
using SimChA.IO;

namespace SimChA.Simulation;

public class MatchSimulator(Random rnd, RefGen refGen, SimParams simParams, FitParams fitParams, EvoParams evoParams)
    : Simulator(rnd, refGen, simParams, fitParams)
{
    private EvoParams EvoParams { get; } = evoParams;

    private double GetExplorationWeight(int tryNo, int totalAttempts, double decay)
    {
        if (totalAttempts < 2 || EvoParams.Decay <= 0 || decay <= 0)
        {
            return 0;
        }

        double progress = (double) tryNo / (totalAttempts - 1);
        double normalizedDecay = Math.Clamp(decay / EvoParams.Decay, 0.0, 1.0);
        return (1.0 - progress) * (1.0 - normalizedDecay);
    }

    private (Karyotype newKar, BaseEventData eventData, int numTries, string signature) GetNewEvent(
        List<CNEventPars> cnEventPars,
        Karyotype currentKar,
        double targetFitness,
        double decay)
    {
        double currentDist = Math.Abs(currentKar.FitnessVal - targetFitness);
        int totalAttempts = EvoParams.MaxTries + 1;

        Karyotype? bestKar = null;
        BaseEventData? bestEventData = null;
        string bestSignature = "";
        double bestScore = double.NegativeInfinity;
        double bestDistImprovement = double.NegativeInfinity;
        bool bestAccepted = false;
        int attemptsPerformed = 0;

        for (int tryNo = 0; tryNo < totalAttempts; tryNo++)
        {
            if (SampleStat.CalcPloidy(currentKar, RefGen) > 32)
            {
                break;
            }

            attemptsPerformed = tryNo + 1;

            var cnEventP = Rnd.PickRndElem(cnEventPars);
            var eventData = Sampling.GenerateCNEventData(Rnd, currentKar, cnEventP);
            if (eventData == null)
            {
                continue;
            }

            var proposedKar = new Karyotype(currentKar);
            eventData.ApplyEvent(proposedKar);
            double proposedFitness = proposedKar.UpdateFitness(RefGen, FitParams);
            double proposedDist = Math.Abs(proposedFitness - targetFitness);

            double distImprovement = currentDist - proposedDist;
            double explorationWeight = GetExplorationWeight(tryNo, totalAttempts, decay);
            double acceptanceProb = Math.Min(1.0, Math.Exp(distImprovement - decay));
            bool accepted = acceptanceProb > Rnd.NextDouble();
            double candidateScore = distImprovement + explorationWeight * acceptanceProb;

            bool tiesBestScore = Math.Abs(candidateScore - bestScore) < double.Epsilon;
            bool tiesBestImprovement = Math.Abs(distImprovement - bestDistImprovement) < double.Epsilon;
            bool preferProposal = candidateScore > bestScore ||
                                 (tiesBestScore && distImprovement > bestDistImprovement) ||
                                 (tiesBestScore && tiesBestImprovement && accepted && !bestAccepted);

            if (preferProposal)
            {
                bestKar = proposedKar;
                bestEventData = eventData;
                bestSignature = cnEventP.Signature;
                bestScore = candidateScore;
                bestDistImprovement = distImprovement;
                bestAccepted = accepted;
            }
        }

        if (bestKar != null && bestEventData != null)
        {
            return (bestKar, bestEventData, Math.Max(attemptsPerformed - 1, 0), bestSignature);
        }

        return (currentKar, CreateSkipEvent(), EvoParams.MaxTries, "");
    }

    protected override (Karyotype childKar, List<CNEventDesc> childEvs) SampleEvents(
        Karyotype currentKar,
        CTreeNode cnChild,
        List<CNEventPars> cnEventPs,
        int mutDepth)
    {
        var childEvs = new List<CNEventDesc>();
        double targetFit = cnChild.Fitness;

        int eventCount = SampleEventCount(cnChild);
        for (int evNo = 1; evNo <= eventCount; evNo++)
        {
            Console.Write($"\rSample {cnChild.CloneId}. Event {evNo}/{eventCount}.".PadRight(80));
            double decay = EvoParams.Decay * evNo / eventCount;
            (var childKar, var eventData, int numTries, string signature) = GetNewEvent(cnEventPs, currentKar, targetFit, decay);
            (string gainedStr, string lostStr) = CalcKaryotypeDiff(currentKar, childKar);
            string karStr = CNEventDesc.PrintKaryotype ? childKar.ToString() : "";
            double newFit = childKar.FitnessVal;
            double dFit = newFit - currentKar.FitnessVal;
            var newEv = new CNEventDesc(eventData, mutDepth + evNo, dFit, newFit, numTries, signature,
                RegionsGained: gainedStr, RegionsLost: lostStr, Karyotype: karStr);
            childEvs.Add(newEv);
            currentKar = childKar;
        }

        return (currentKar, childEvs);
    }
}
