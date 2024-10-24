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

    private double CalculateAcceptance(double newFitness, double oldFitness, double temperature, double mutationRate = 1)
    {
        var mutPart = Math.Log(mutationRate);
        if (EvoParams.WithFitness)
        {
            return Math.Min(0, mutPart);
        }
        var fitPart = (newFitness - oldFitness)/Math.Abs(oldFitness);
        if (EvoParams.SimulatedAnnealing)
        {
            fitPart /= temperature;
        }
        else
        {
            fitPart *= EvoParams.ThetaFitness;
        }
        return Math.Min(0, fitPart + mutPart);
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
        var currentTemp = EvoParams.Temperature;

        for (int i = 0; i < EvoParams.NumIterations; i++)
        {
            Console.Write($"\rSample {sample.SampleId}. Iteration {i+1}/{EvoParams.NumIterations}; Event Count {currentEvents.Count}.".PadRight(80));
            // Generate a new event and correspondingly add to list
            var newEvent = GetNewEvent(sample, new Karyotype(kar));
            var proposedFitness = GetFitness(new Karyotype(kar), newEvent);
            var acceptProb = CalculateAcceptance(proposedFitness, currentFitness, EvoParams.MutationRate, currentTemp);
            if (acceptProb >= Math.Log(Rnd.NextDouble()))
            {
                currentFitness = proposedFitness;
                currentEvents.Add(newEvent);
                // Apply the new event to the clone
                newEvent.ApplyEvent(kar);
                kar.UpdateFitness(GenRef, FitnessParams);
            }
            currentTemp *= EvoParams.CoolingRate;
        }
        return currentEvents;
    }

    private List<BaseEventData> EvolveInEvents(Sample sample, Karyotype kar, int mutCount)
    {
        var currentEvents = new List<BaseEventData>();
        var currentFitness = Fitness.Calculate(new Karyotype(kar), GenRef, FitnessParams);
        var currentTemp = EvoParams.Temperature;
        //var bestFitness = currentFitness;

        for (int i = 0; i < EvoParams.NumIterations && currentEvents.Count < mutCount; i++)
        {
            Console.Write($"\rSample {sample.SampleId}. Iteration {i+1}/{EvoParams.NumIterations}; Event Count {currentEvents.Count}/{mutCount}.".PadRight(80));
            // Generate a new event and correspondingly add to list
            var newEvent = GetNewEvent(sample, new Karyotype(kar));
            var proposedFitness = GetFitness(new Karyotype(kar), newEvent);
            var acceptProb = CalculateAcceptance(proposedFitness, currentFitness, currentTemp);
            if (acceptProb >= Math.Log(Rnd.NextDouble()))
            {
                currentFitness = proposedFitness;
                currentEvents.Add(newEvent);
                // Apply the new event to the clone
                newEvent.ApplyEvent(kar);
                kar.UpdateFitness(GenRef, FitnessParams);
            }
            currentTemp *= EvoParams.CoolingRate;
        }
        return currentEvents;
    }

    private void ApplyEvolutionRec(Sample sample, CloneIn node, IReadOnlyDictionary<int, 
        List<CloneIn>> clones, int eventCount)
    {
        foreach (var child in clones[node.CloneId])
        {
            var childKar = new Karyotype(sample.Kars[node.CloneId]);
            // copy of karyotype for printing out the events & their individual effects
            var dummyKar = new Karyotype(sample.Kars[node.CloneId]);
            sample.Kars[child.CloneId] = childKar;
            var childEvs = new List<CNEventDesc>();
            sample.EventDescs[child.CloneId] = childEvs;
            
            double oldFitness = Fitness.Calculate(childKar, GenRef, FitnessParams);

            var bestEvents = EvoParams.EvolveInTime 
                ? EvolveInTime(sample, childKar)
                : EvolveInEvents(sample, childKar, child.Distance);
            Console.WriteLine("Fetching the sampled events and calculating fitness changes");
            
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
