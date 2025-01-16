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
        ApplyEvolutionRec(sample, root, childLookUp, 0);
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

    private (List<CNEventDesc>, Karyotype) EvolveInEvents(Sample sample, Karyotype kar, int eventCount, int mutCount)
    {
        var childEvs = new List<CNEventDesc>();
        double currentFitness = Fitness.Calculate(kar, GenRef, FitnessParams);
        for (int j = 0; 
             j < mutCount * EvoParams.MaxTries && childEvs.Count < mutCount;
             j = Math.Max(j + 1, childEvs.Count * EvoParams.MaxTries))
        {
            var proposedKar = new Karyotype(kar);
            Console.Write($"\rSample {sample.SampleId}. Iteration {j/mutCount}/{mutCount};".PadRight(80));
            var cnEventPars = EvoParams.EventCost 
                ? GetModifiedEventPars(sample.EventPars) 
                : sample.EventPars;
            var newEvent = GetNewEvent(cnEventPars, proposedKar);
            newEvent.ApplyEvent(proposedKar);
            double proposedFitness = proposedKar.UpdateFitness(GenRef, FitnessParams);
            double dFit = proposedFitness - currentFitness;
            if (EvoParams.WithFitness && dFit < Math.Log(Rnd.NextDouble()))
            {
                continue;
            }
            var abberation = new CNEventDesc(
                newEvent.EventType,
                eventCount + childEvs.Count,
                newEvent.ToString(),
                dFit,
                proposedFitness,
                childEvs.Count + 1);
            childEvs.Add(abberation);
            kar = proposedKar;
            currentFitness = proposedFitness;
        }
        if (childEvs.Count < mutCount)
        {
            throw new Exception("Failed to generate the required number of events for sample " + sample.SampleId);
        }

        return (childEvs, kar);
    }
    
    private void ApplyEvolutionRec(
        Sample sample,
        CloneIn node,
        IReadOnlyDictionary<string, List<CloneIn>> clones,
        int eventCount)
    {
        foreach (var child in clones[node.CloneId])
        {
            // copy of karyotype for printing out the events & their individual effects
            // Start with the tetraploid state
            if (EvoParams.TetraploidStart)
            {
                sample.Kars[node.CloneId].ApplyWGD();
            }

            var (newEvents, newKaryotype) = EvolveInEvents(sample, sample.Kars[node.CloneId], eventCount, child.Distance);
            sample.EventDescs[node.CloneId] = newEvents;
            sample.Kars[node.CloneId] = newKaryotype;
            
            // The sample's clone should have its distance updated
            int cloneIndex = sample.Clones.FindIndex(c => c.CloneId == child.CloneId);
            if (cloneIndex != -1)
            {
                var updatedClone = sample.Clones[cloneIndex] with { Distance = eventCount + sample.EventDescs[node.CloneId].Count };
                sample.Clones[cloneIndex] = updatedClone;
            }
            else
            {
                throw new Exception("Error in Evolver. ApplyEvolutionRec: Clone not found in sample.");
            }

            if (child.CloneId != node.CloneId)
            {
                ApplyEvolutionRec(sample, child, clones, eventCount + child.Distance);
            }
        }
    }
}