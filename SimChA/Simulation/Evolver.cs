using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;
using SimChA.IO;

namespace SimChA.Simulation;

public class Evolver
{
    protected readonly EvoParams EvoParams;
    protected readonly FitnessParams FitnessParams;
    protected readonly GenRef GenRef;
    protected readonly Random Rnd;
    protected List<double>? EventTimes;
    protected List<CNEventPars>? PostWGDEventPars;
    protected List<CNEventPars>? PreWGDEventPars;

    public Evolver(
        Random rnd,
        GenRef genRef,
        FitnessParams fitnessParams,
        EvoParams evoParams)
    {
        Rnd = rnd;
        GenRef = genRef;
        FitnessParams = fitnessParams;
        EvoParams = evoParams;
    }

    public void EvolveSample(Sample sample)
    {
        if (sample.EventPars == null || sample.EventPars.Count == 0)
        {
            throw new Exception("No events to sample from.");
        }

        if (sample.Signatures.Count == 2
            && sample.Signatures.TryGetValue("PostWGD", out var postWGDSig)
            && sample.Signatures.TryGetValue("PreWGD", out var preWGDSig))
        {
            PreWGDEventPars = Converters.NormalizeEvents(preWGDSig.Events);
            PostWGDEventPars = Converters.NormalizeEvents(postWGDSig.Events);
            if (EvoParams.TetraploidStart)
            {
                PreWGDEventPars = PostWGDEventPars;
            }
        }

        var (root, childLookUp) = CloneComp.CreateLookUp(sample.Clones);
        sample.Kars[root.CloneId] = new Karyotype(GenRef, sample.Sex);
        ApplyEvolutionRec(sample, root, childLookUp, 1);
    }

    private bool DidMutate(double mean)
        => Rnd.NextDouble() < 1 - Math.Exp(-mean);

    private double CalculateLogAcceptance(double newFitness, double oldFitness)
    {
        if (!EvoParams.WithFitness)
        {
            return 0;
        }
        var fitPart = EvoParams.ThetaFitness * (newFitness - oldFitness);
        return Math.Min(0, fitPart);
    }

    private double GetFitness(Karyotype kar, List<BaseEventData> eventData)
    {
        foreach (var ev in eventData)
        {
            ev.ApplyEvent(kar);
        }
        return kar.UpdateFitness(GenRef, FitnessParams);
    }

    private int GetEventCount(Karyotype kar, bool hasDoubled = false)
    {
        int nEvents;
        /*var mu = EvoParams.DynamicMutRate
                ? EvoParams.MutationRate * CNProfile.CalcPloidy(kar, GenRef) / 2.0
                : EvoParams.MutationRate; */
        var mu = EvoParams.DynamicMutRate && hasDoubled
            ? EvoParams.MutationRate * EvoParams.WGDAccelerationFactor
            : EvoParams.MutationRate;
        if (EvoParams.EventBlock)
        {
            if (EvoParams.StepDistribution != Distribution.Poisson)
            {
                throw new Exception("Invalid distribution for event block.");
            }

            nEvents = Sampling.SampleDistInt(Rnd, EvoParams.StepDistribution, mu);
        }
        else
        {
            nEvents = DidMutate(mu) ? 1 : 0;
        }

        return nEvents;
    }

    private List<CNEventPars> GetModifiedEventPars(List<CNEventPars> pars)
    {
        var newPars = new List<CNEventPars>(pars);
        //var factor = CNProfile.CalcPloidy(kar, GenRef)/2.0;
        foreach (var e in pars)
        {
            var newProb = e.Type switch
            {
                CNEventType.ChromDeletion => Math.Max(0, e.Prob * 4), //* EvoParams.ChromLossEnhancementFactor),
                CNEventType.ArmDeletion => Math.Max(0, e.Prob * 2), //* EvoParams.ChromLossEnhancementFactor),
                _ => e.Prob
            };
            newPars[newPars.IndexOf(e)] = e with { Prob = newProb };
        }

        var normalized = Converters.NormalizeEvents(newPars);
        return normalized;
    }

    private List<BaseEventData> GetNewEvents(
        List<CNEventPars> eventPars,
        Karyotype kar,
        int nEvents,
        bool wgdPositive = false)
    {
        var sampledEvents = new List<BaseEventData>();
        var iTries = 0;
        var pars = EvoParams.EventCost && wgdPositive
            ? GetModifiedEventPars(eventPars)
            : eventPars;
        for (var i = 0; i < nEvents && iTries < EvoParams.MaxTries;)
        {
            var cnEventP = Rnd.PickRndElem(pars);
            var eventData = Sampling.GenerateCNEventData(Rnd, kar, cnEventP);

            if (eventData != null)
            {
                sampledEvents.Add(eventData);
                i++;
            }
            else
            {
                iTries++;
            }
        }

        if (iTries >= EvoParams.MaxTries)
        {
            throw new Exception("MaxTries exceeded in Evolver.cs. Could not generate new events.");
        }

        return sampledEvents;
    }

    private List<BaseEventData> EvolveInContinuousTime(Sample sample, Karyotype kar)
    {
        var currentEvents = new List<BaseEventData>();
        var currentFitness = Fitness.Calculate(new Karyotype(kar), GenRef, FitnessParams);
        var timeList = new List<double> { 0.0 };
        EventTimes = new List<double>();
        var hasDoubled = false;
        var eventPars = PreWGDEventPars ?? sample.EventPars;
        var wgdCount = 0;
        while (timeList.Last() < EvoParams.MaxTime)
        {
            // Sample the new time for the event
            var u = Rnd.NextDouble();
            var mu = EvoParams.DynamicMutRate && hasDoubled
                ? EvoParams.MutationRate * EvoParams.WGDAccelerationFactor
                : EvoParams.MutationRate;
            var tNew = timeList.Last() - Math.Log(u) / mu;
            if (tNew > EvoParams.MaxTime)
            {
                break;
            }

            timeList.Add(tNew);
            // Generate a new event and correspondingly add to list
            var newEvents = GetNewEvents(eventPars, new Karyotype(kar), 1, hasDoubled);
            if (newEvents.Count > 1)
            {
                throw new Exception("Continuous time evolution should only sample one event at a time.");
            }

            if (newEvents.Count == 0)
            {
                continue;
            }

            var ev = newEvents[0];
            var proposedFitness = GetFitness(new Karyotype(kar), newEvents);
            var acceptProb = CalculateLogAcceptance(proposedFitness, currentFitness);
            if (acceptProb >= Math.Log(Rnd.NextDouble()))
            {
                currentFitness = proposedFitness;
                currentEvents.Add(ev);
                ev.ApplyEvent(kar);
                EventTimes.Add(tNew);
                kar.UpdateFitness(GenRef, FitnessParams);
                if (!hasDoubled && ev.EventType == CNEventType.WholeGenomeDoubling)
                {
                    eventPars = PostWGDEventPars ?? sample.EventPars;
                    hasDoubled = true;
                }

                if (ev.EventType == CNEventType.WholeGenomeDoubling)
                {
                    wgdCount += 1;
                }
            }

            // Short-circuit the evolution process if we already have too many WGDs
            if (wgdCount >= 6)
            {
                break;
            }
        }

        return currentEvents;
    }

    private List<BaseEventData> EvolveInTime(Sample sample, Karyotype kar)
    {
        var currentEvents = new List<BaseEventData>();
        var currentFitness = Fitness.Calculate(new Karyotype(kar), GenRef, FitnessParams);
        EventTimes = new List<double>();
        var eventPars = sample.EventPars;
        for (var i = 0; i < EvoParams.MaxTime; i++)
        {
            Console.Write(
                $"\rSample {sample.SampleId}. Iteration {i + 1}/{EvoParams.MaxTime}; Event Count {currentEvents.Count}."
                    .PadRight(80));
            // Generate a new event and correspondingly add to list
            // Want to sample a number of events.
            var nEvents = GetEventCount(kar);
            var newEvents = GetNewEvents(eventPars, new Karyotype(kar), nEvents);
            if (newEvents.Count == 0)
            {
                continue;
            }

            var proposedFitness = GetFitness(new Karyotype(kar), newEvents);
            var acceptProb = CalculateLogAcceptance(proposedFitness, currentFitness);
            if (acceptProb >= Math.Log(Rnd.NextDouble()))
            {
                currentFitness = proposedFitness;
                foreach (var ev in newEvents)
                {
                    currentEvents.Add(ev);
                    ev.ApplyEvent(kar);
                    EventTimes.Add(i);
                }

                kar.UpdateFitness(GenRef, FitnessParams);
            }
        }

        return currentEvents;
    }

    private List<BaseEventData> EvolveInEvents(Sample sample, Karyotype kar, int mutCount)
    {
        var currentEvents = new List<BaseEventData>();
        var currentFitness = Fitness.Calculate(new Karyotype(kar), GenRef, FitnessParams);

        var nSteps = mutCount; //GetNumSteps(mutCount, kar);
        var eventPars = sample.EventPars;
        var i = 0;
        var hasDoubled = false;
        for (; i < nSteps;)
        {
            Console.Write($"\rSample {sample.SampleId}. Iteration {i + 1}/{nSteps};".PadRight(80));
            // Generate a new event and correspondingly add to list
            var nEvents = 1; // GetEventCount(kar);
            var newEvents = GetNewEvents(eventPars, new Karyotype(kar), nEvents);
            if (newEvents.Count != 1)
            {
                continue;
            }

            var proposedFitness = GetFitness(new Karyotype(kar), newEvents);
            var acceptProb = CalculateLogAcceptance(proposedFitness, currentFitness);
            if (acceptProb >= Math.Log(Rnd.NextDouble()))
            {
                currentFitness = proposedFitness;
                foreach (var ev in newEvents)
                {
                    currentEvents.Add(ev);
                    ev.ApplyEvent(kar);
                    if (!hasDoubled && ev.EventType == CNEventType.WholeGenomeDoubling)
                    {
                        hasDoubled = true;
                    }
                }

                kar.UpdateFitness(GenRef, FitnessParams);
                i++;
            }
        }

        return currentEvents;
    }

    private List<BaseEventData> GetEvents(Sample sample, CloneIn clone)
    {
        var childKar = sample.Kars[clone.CloneId];
        if (!EvoParams.EvolveInTime)
        {
            return EvolveInEvents(sample, childKar, clone.Distance);
        }
        if (!EvoParams.ContinuousTime)
        {
            return EvolveInTime(sample, childKar);
        }
        return EvolveInContinuousTime(sample, childKar);
    }

    private void ApplyEvolutionRec(
        Sample sample,
        CloneIn node,
        IReadOnlyDictionary<string, List<CloneIn>> clones,
        int eventCount)
    {
        foreach (var child in clones[node.CloneId])
        {
            var childKar = new Karyotype(sample.Kars[node.CloneId]);
            // copy of karyotype for printing out the events & their individual effects
            var dummyKar = new Karyotype(sample.Kars[node.CloneId]);
            // Start with the tetraploid state
            if (EvoParams.TetraploidStart)
            {
                childKar.ApplyWGD();
                dummyKar.ApplyWGD();
            }

            sample.Kars[child.CloneId] = childKar;
            var childEvs = new List<CNEventDesc>();
            sample.EventDescs[child.CloneId] = childEvs;
            var oldFitness = Fitness.Calculate(childKar, GenRef, FitnessParams);
            var bestEvents = GetEvents(sample, child);

            for (var mutNo = 0; mutNo < bestEvents.Count; mutNo++)
            {
                Console.Write($"\rSample {sample.SampleId}." +
                              $" Clone {child.CloneId}" +
                              $" Event {mutNo + 1}/{bestEvents.Count}.");
                var eventData = bestEvents[mutNo];
                eventData.ApplyEvent(dummyKar);
                var newFitness = dummyKar.UpdateFitness(GenRef, FitnessParams);
                var dFit = newFitness - oldFitness;
                var time = EvoParams.EvolveInTime && EventTimes != null ? EventTimes[mutNo] : 0;
                var abberation = new CNEventDesc(
                    eventData.EventType,
                    eventCount + mutNo,
                    eventData.ToString(),
                    dFit,
                    newFitness,
                    time);
                childEvs.Add(abberation);
                oldFitness = newFitness;
            }

            // The sample's clone should have its distance updated
            var cloneIndex = sample.Clones.FindIndex(c => c.CloneId == child.CloneId);
            if (cloneIndex != -1)
            {
                var updatedClone = sample.Clones[cloneIndex] with { Distance = eventCount + bestEvents.Count };
                sample.Clones[cloneIndex] = updatedClone;
            }
            else
            {
                throw new Exception("Error in Evolver.ApplyEvolutionRec: Clone not found in sample.");
            }
            if (child.CloneId != node.CloneId)
            {
                ApplyEvolutionRec(sample, child, clones, eventCount + child.Distance);
            }
        }
    }
}