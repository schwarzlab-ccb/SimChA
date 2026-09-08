using SimChA.Computation;
using SimChA.Data;
using SimChA.EventData;
using SimChA.IO;

namespace SimChA.Simulation;

public class EvoSimulator(Random rnd, RefGen refGen, SimParams simParams, FitParams fitParams, EvoParams evoParams)
    : Simulator(rnd, refGen, simParams, fitParams)
{
    private EvoParams EvoParams { get; } = evoParams;

    /// <summary>How a mutation slot ended. Recorded per event; see SlotResult.</summary>
    private const string OutcomeAccepted = "Accepted";
    private const string OutcomeImpossible = "Impossible";
    private const string OutcomeRejected = "Rejected";
    private const string OutcomePloidyCeiling = "PloidyCeiling";

    /// <summary>Ploidy past which no further event is attempted and the slot is skipped.</summary>
    private const double MaxPloidy = 32;

    /// <summary>
    /// One mutation slot's result, with the tries it spent broken down by what consumed them.
    ///
    /// The breakdown is instrumentation, not behaviour: impossibility and inviability both still
    /// `continue`, exactly as before. It exists because their rates are otherwise unobservable --
    /// the retry loop redraws the event type, so a type that is impossible on this karyotype is
    /// silently resampled away rather than consuming the slot the way the base Simulator's single
    /// draw does. The impossibility rate is what a slot-consuming exit would turn into skips, and
    /// that number has to be measured before the exit can be changed, because it also sets how far
    /// the realized event count would fall below the configured RateMean.
    /// </summary>
    private readonly record struct SlotResult(
        Karyotype NewKar,
        BaseEventData EventData,
        int NumRejections,
        int NumImpossible,
        int NumInviable,
        string Outcome,
        string AttemptedType,
        string Signature);

    private SlotResult GetNewEvent(List<CNEventPars> cnEventPars, Karyotype currentKar)
    {
        int impossible = 0, inviable = 0, rejected = 0;
        string attemptedType = "";
        for (int tryNo = 0; tryNo <= EvoParams.MaxTries; tryNo++)
        {
            if (SampleStat.CalcPloidy(currentKar, RefGen) > MaxPloidy)
            {
                return new SlotResult(currentKar, CreateSkipEvent(), rejected, impossible, inviable,
                                      OutcomePloidyCeiling, attemptedType, "");
            }
            var cnEventP = Rnd.PickRndElem(cnEventPars);
            attemptedType = cnEventP.Type.ToString();
            var eventData = Sampling.GenerateCNEventData(Rnd, currentKar, cnEventP);
            if (eventData == null)
            {
                impossible++;
                continue;
            }

            var proposedKar = new Karyotype(currentKar);
            eventData.ApplyEvent(proposedKar);
            // Before the fitness is even computed: an essential gene at zero copies is inviable, so
            // the proposal is rejected outright rather than priced. See Fitness.AnyEssentialLost.
            if (FitParams.ProhibitEssentialLoss &&
                Fitness.AnyEssentialLost(RefGen.SexGeneLists[(int) proposedKar.Sex][(int) GeneLT.Ess],
                                         proposedKar.GeneCounts[(int) GeneLT.Ess]))
            {
                inviable++;
                continue;
            }
            double proposedFitness = proposedKar.UpdateFitness(RefGen, FitParams);
            if (Fitness.AcceptProb(proposedFitness - currentKar.FitnessVal, EvoParams.Acceptance) > Rnd.NextDouble())
            {
                return new SlotResult(proposedKar, eventData, rejected, impossible, inviable,
                                      OutcomeAccepted, attemptedType, cnEventP.Signature);
            }
            rejected++;
        }
        return new SlotResult(currentKar, CreateSkipEvent(), rejected, impossible, inviable,
                              OutcomeRejected, attemptedType, "");
    }
    
    protected override (Karyotype childKar, List<CNEventDesc> childEvs) SampleEvents(
        Karyotype currentKar,
        CTreeNode cnChild,
        List<CNEventPars> cnEventPs,
        int mutDepth,
        int eventCount)
    {
        var childEvs = new List<CNEventDesc>();

        for (int evNo = 1; evNo <= eventCount; evNo++)
        {
            Console.Write($"\rSample {cnChild.CloneId}. Event {evNo}/{eventCount}.".PadRight(80));
            var slot = GetNewEvent(cnEventPs, currentKar);
            var childKar = slot.NewKar;
            (string gainedStr, string lostStr) = CalcKaryotypeDiff(currentKar, childKar);
            string karStr =  CNEventDesc.PrintKaryotype ? childKar.ToString() : "";
            double newFit = childKar.FitnessVal;
            double dFit = newFit - currentKar.FitnessVal;
            var newEv = new CNEventDesc(slot.EventData, mutDepth + evNo, dFit, newFit,
                slot.NumRejections, slot.NumImpossible, slot.NumInviable, slot.Outcome,
                slot.AttemptedType, slot.Signature,
                RegionsGained: gainedStr, RegionsLost: lostStr, Karyotype: karStr);
            childEvs.Add(newEv);
            currentKar = childKar;
        }
        return (currentKar, childEvs);
    }
}
