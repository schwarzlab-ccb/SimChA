using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;
using SimChA.IO;
using EDists = Extreme.Statistics.Distributions;

namespace SimChA.Simulation;

public class Evolver
{
    protected readonly FitnessParams FitnessParams;
    protected readonly EvoParams EvoParams;
    protected readonly FileIO FileIO;
    protected readonly Random Rnd;
    protected readonly GenRef GenRef;
    protected int Counter;
    protected List<CNEventPars>? EventPars = null;

    public Evolver(
        Random rnd,
        GenRef genRef,
        FitnessParams fitnessParams, 
        EvoParams evoParams,
        FileIO fileIO)
    {
        Rnd = rnd;
        GenRef = genRef;
        FitnessParams = fitnessParams;
        EvoParams = evoParams;
        FileIO = fileIO;
    }

    public void EvolveSample(Sample sample)
    {
        if (sample.EventPars == null || !sample.EventPars.Any())
        {
            throw new Exception("No events to sample from.");
        }
        Counter = 1;
        EventPars = sample.EventPars;
        var (root, childLoopUp) = CloneComp.CreateLookUp(sample.Clones);
        sample.Kars[root.CloneId] = new Karyotype(GenRef, sample.Sex);
        ApplyEvolutionRec(sample, root, childLoopUp, 1);
    }

    private double CalculatePotential(double proposedFitness)
    {
        double fitnessPotential = EvoParams.ThetaFitness * proposedFitness;
        return fitnessPotential;
    }    

    public double GetFitness(Karyotype kar, BaseEventData eventData)
    {
        eventData.ApplyEvent(kar);
        return kar.UpdateFitness(GenRef, FitnessParams);
    }

    private BaseEventData GetNewEvent(Sample sample, Karyotype kar)
    {
        var cnEventP = Rnd.PickRndElem(sample.EventPars);
        var eventData = Sampling.GenerateCNEventData(Rnd, kar, cnEventP);
        return eventData ?? throw new Exception("Failed to generate new event data.");
    }

    private List<BaseEventData> EvolveInTime(Sample sample, Karyotype kar)
    {
        var currentEvents = new List<BaseEventData>();
        var currentFitness = Fitness.Calculate(new Karyotype(kar), GenRef, FitnessParams);
        var currentPotential = CalculatePotential(currentFitness);

        var bestFitness = currentFitness;

        for (int i = 0; i < EvoParams.NumIterations; i++)
        {
            Console.Write($"\rSample {sample.SampleId}. Event {currentEvents.Count+1}/{EvoParams.NumIterations}.".PadRight(80));
            // Generate a new event and correspondingly add to list
            var newEvent = GetNewEvent(sample, new Karyotype(kar));
            var proposedFitness = GetFitness(new Karyotype(kar), newEvent);
            var proposedPotential = CalculatePotential(proposedFitness);
            var acceptProb = Math.Min(0, proposedPotential - currentPotential) * EvoParams.MutationRate;
            if (acceptProb >= Math.Log(Rnd.NextDouble()))
            {
                currentPotential = proposedPotential;
                currentEvents.Add(newEvent);
                if (proposedFitness > bestFitness)
                {
                    bestFitness = proposedFitness;
                }
                // Apply the new event to the clone
                newEvent.ApplyEvent(kar);
            }
        }
        return currentEvents;
    }

    private List<BaseEventData> EvolveInEvents(Sample sample, Karyotype kar, int mutCount)
    {
        var currentEvents = new List<BaseEventData>();
        var currentFitness = Fitness.Calculate(new Karyotype(kar), GenRef, FitnessParams);
        var currentPotential = CalculatePotential(currentFitness);

        var bestFitness = currentFitness;

        for (int i = 0; i < EvoParams.NumIterations && currentEvents.Count < mutCount; i++)
        {
            Console.Write($"\rSample {sample.SampleId}. Event {currentEvents.Count+1}/{EvoParams.NumIterations}.".PadRight(80));
            // Generate a new event and correspondingly add to list
            var newEvent = GetNewEvent(sample, new Karyotype(kar));
            var proposedFitness = GetFitness(new Karyotype(kar), newEvent);
            var proposedPotential = CalculatePotential(proposedFitness);
            var acceptProb = Math.Min(0, proposedPotential - currentPotential);
            if (acceptProb >= Math.Log(Rnd.NextDouble()))
            {
                currentPotential = proposedPotential;
                currentEvents.Add(newEvent);
                if (proposedFitness > bestFitness)
                {
                    bestFitness = proposedFitness;
                }
                // Apply the new event to the clone
                newEvent.ApplyEvent(kar);
            }
        }
        return currentEvents;
    }

    private void ApplyEvolutionRec(Sample sample, CloneIn node, IReadOnlyDictionary<int, 
        List<CloneIn>> clones, int eventCount)
    {
        foreach (var child in clones[node.CloneId])
        {
            var childKar = new Karyotype(sample.Kars[node.CloneId]);
            sample.Kars[child.CloneId] = childKar;
            var childEvs = new List<CNEventDesc>();
            sample.EventDescs[child.CloneId] = childEvs;
            
            double oldFitness = Fitness.Calculate(childKar, GenRef, FitnessParams);

            var bestEvents = EvoParams.EvolveInTime 
                ? EvolveInTime(sample, childKar)
                : EvolveInEvents(sample, childKar, child.Distance);
            Console.WriteLine("Fetching the sampled events and calculating fitness changes");
            var dummyKar = new Karyotype(sample.Kars[node.CloneId]);
            for (int mutNo = 0; mutNo < bestEvents.Count; mutNo++)
            {
                Console.Write($"\rSample {sample.SampleId}. Clone {Counter}/{clones.Count}. Event {mutNo + 1}/{bestEvents.Count}.");
                var eventData = bestEvents[mutNo];
                eventData.ApplyEvent(dummyKar);
                double newFitness = dummyKar.UpdateFitness(GenRef, FitnessParams);
                double dFit = newFitness - oldFitness;
                var abberation = new CNEventDesc(eventData.EventType, eventCount + mutNo, eventData.ToString(), dFit, newFitness);
                childEvs.Add(abberation);
                oldFitness = newFitness;
            }

            Counter++;
            if (child.CloneId != node.CloneId)
            {
                ApplyEvolutionRec(sample, child, clones, eventCount + child.Distance);
            }
        }
    }
}
