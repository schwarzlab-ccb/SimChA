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

        var (root, childLookUp) = CloneComp.CreateLookUp(sample.Clones);
        sample.Kars[root.CloneId] = new Karyotype(GenRef, sample.Sex);
        ApplyEvolutionRec(sample, root, childLookUp, 1);
    }

    private double CalcFitnessChange(double newFitness, double oldFitness)
        => EvoParams.WithFitness ? newFitness - oldFitness : 0;
    
    private double ApplyEvents(Karyotype kar, List<BaseEventData> eventData)
    {
        foreach (var ev in eventData)
        {
            ev.ApplyEvent(kar);
        }
        return kar.UpdateFitness(GenRef, FitnessParams);
    }
    
    private List<CNEventPars> GetModifiedEventPars(List<CNEventPars> pars)
    {
        var newPars = new List<CNEventPars>(pars);
        //var factor = CNProfile.CalcPloidy(kar, GenRef)/2.0;
        foreach (var e in pars)
        {
            double newProb = e.Type switch
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

    private BaseEventData GetNewEvent(List<CNEventPars> cnEventPars, Karyotype kar)
    {
        for (int tryNo = 0; tryNo <= EvoParams.MaxTries; tryNo++)
        {
            var cnEventP = Rnd.PickRndElem(cnEventPars);
            var eventData = Sampling.GenerateCNEventData(Rnd, kar, cnEventP);
            if (eventData != null)
            {
                return eventData;
            }
        }
        throw new Exception($"Could not generate a new event. In {EvoParams.MaxTries} tries.");
    }

    private List<BaseEventData> GetNewEvents(
        List<CNEventPars> eventPars,
        Karyotype kar,
        int nEvents,
        bool wgdPositive = false)
    {
        var sampledEvents = new List<BaseEventData>();
        var cnEventPars = EvoParams.EventCost && wgdPositive ? GetModifiedEventPars(eventPars) : eventPars;
        for (int i = 0; i < nEvents; i++)
        {
            sampledEvents.Add(GetNewEvent(cnEventPars, kar));
        }
        return sampledEvents;
    }

    private List<BaseEventData> EvolveInEvents(Sample sample, Karyotype kar, int mutCount)
    {
        var currentEvents = new List<BaseEventData>();
        double currentFitness = Fitness.Calculate(new Karyotype(kar), GenRef, FitnessParams);
        int nSteps = mutCount; //GetNumSteps(mutCount, kar);
        var eventPars = sample.EventPars;
        int i = 1;
        bool hasDoubled = false;
        while (i < nSteps)
        {
            Console.Write($"\rSample {sample.SampleId}. Iteration {i + 1}/{nSteps};".PadRight(80));
            // Generate a new event and correspondingly add to list
            int nEvents = 1; // GetEventCount(kar);
            var newEvents = GetNewEvents(eventPars, new Karyotype(kar), nEvents);
            if (newEvents.Count != 1)
            {
                continue;
            }

            double proposedFitness = ApplyEvents(new Karyotype(kar), newEvents);
            double acceptProb = CalcFitnessChange(proposedFitness, currentFitness);
            if (acceptProb < Math.Log(Rnd.NextDouble()))
            {
                continue;
            }
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
        return currentEvents;
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
            double oldFitness = Fitness.Calculate(childKar, GenRef, FitnessParams);
            var bestEvents = EvolveInEvents(sample, childKar, child.Distance);
            
            for (int mutNo = 0; mutNo < bestEvents.Count; mutNo++)
            {
                string progressStr = $"\rSample {sample.SampleId}." +
                                     $" Clone {child.CloneId}." +
                                     $" Event {mutNo + 1}/{bestEvents.Count}.";
                Console.Write(progressStr.PadRight(80));
                var eventData = bestEvents[mutNo];
                eventData.ApplyEvent(dummyKar);
                double newFitness = dummyKar.UpdateFitness(GenRef, FitnessParams);
                double dFit = newFitness - oldFitness;
                var abberation = new CNEventDesc(
                    eventData.EventType,
                    eventCount + mutNo,
                    eventData.ToString(),
                    dFit,
                    newFitness,
                    mutNo + 1);
                childEvs.Add(abberation);
                oldFitness = newFitness;
            }

            // The sample's clone should have its distance updated
            int cloneIndex = sample.Clones.FindIndex(c => c.CloneId == child.CloneId);
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