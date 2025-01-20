using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;
using SimChA.IO;

namespace SimChA.Simulation;

public class Evolver
{
    private EvoParams EvoParams { get; }
    private FitnessParams FitnessParams { get; }
    private GenRef GenRef { get; }
    private Random Rnd { get; }

    public Evolver(Random rnd, GenRef genRef, FitnessParams fitnessParams, EvoParams evoParams)
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
        // Start with the tetraploid state
        if (EvoParams.TetraploidStart)
        {
            sample.Kars[root.CloneId].ApplyWGD();
            sample.Kars[root.CloneId].UpdateFitness(GenRef, FitnessParams);
        }
        ApplyEvolutionRec(sample, root, childLookUp, 0);
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
    
    private void ApplyEvolutionRec(
        Sample sample,
        CloneIn node,
        IReadOnlyDictionary<string, List<CloneIn>> clones,
        int eventCount)
    {
        foreach (var child in clones[node.CloneId])
        {
            int mutCount = child.Distance;
            var childEvs = new List<CNEventDesc>();
            double currentFitness = Fitness.Calculate(sample.Kars[node.CloneId], GenRef, FitnessParams);
            bool hasWGD = false;
            for (int j = 0; 
                 j < mutCount * EvoParams.MaxTries && childEvs.Count < mutCount;
                 j = Math.Max(j + 1, childEvs.Count * EvoParams.MaxTries))
            {
                Console.Write($"\rSample {sample.SampleId}. Iteration {j}/{mutCount};".PadRight(80));
                
                var proposedKar = new Karyotype(sample.Kars[node.CloneId]);
                var cnEventPars = GetEventPars(sample.EventPars, hasWGD);
                var newEvent = GetNewEvent(cnEventPars, proposedKar);
                newEvent.ApplyEvent(proposedKar);
                double proposedFitness = proposedKar.UpdateFitness(GenRef, FitnessParams);
                double dFit = proposedFitness - currentFitness;
                if (EvoParams.WithFitness && Math.Exp(dFit) < Rnd.NextDouble())
                {
                    continue;
                }
                var abberation = new CNEventDesc(
                    newEvent.EventType,
                    eventCount + childEvs.Count + 1,
                    newEvent.ToString(),
                    dFit,
                    proposedFitness,
                    childEvs.Count + 1);
                hasWGD |= newEvent.EventType == CNEventType.WholeGenomeDoubling;
                childEvs.Add(abberation);
                sample.Kars[node.CloneId] = proposedKar;
                currentFitness = proposedFitness;
		if (proposedKar.GenomeLen() == 0)
		{
			break;
		}
            }
            if (childEvs.Count < mutCount && sample.Kars[node.CloneId].GenomeLen() > 0)
            {
                throw new Exception("Failed to generate the required number of events for sample " + sample.SampleId);
            }
            sample.EventDescs[node.CloneId] = childEvs;
            if (child.CloneId != node.CloneId)
            {
                ApplyEvolutionRec(sample, child, clones, eventCount + child.Distance);
            }
        }
    }
}
