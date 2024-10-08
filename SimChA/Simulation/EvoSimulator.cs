using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;
using SimChA.IO;
using EDists = Extreme.Statistics.Distributions;

namespace SimChA.Simulation;

public class EvoSimulator : Simulator
{
    private FitnessParams FitnessParams { get; }
    private MCParams McParams { get; }
    private FileIO FileIO { get; }
    public EvoSimulator(
        Random rnd,
        GenRef genRef,
        FitnessParams fitnessParams, 
        MCParams mCParams,
        FileIO fileIO) : base(rnd, genRef)
    {
        FitnessParams = fitnessParams;
        McParams = mCParams;
        FileIO = fileIO;
    }

    public override void EvolveSample(Sample sample)
    {
        if (sample.EventPars == null || !sample.EventPars.Any())
        {
            throw new Exception("No events to sample from.");
        }
        Counter = 1;
        var (root, childLoopUp) = CloneComp.CreateLookUp(sample.Clones);
        sample.Kars[root.CloneId] = new Karyotype(GenRef, sample.Sex);
        ApplyEvolutionRec(sample, root, childLoopUp, 1);
    }

    public double GetEventPotential(List<BaseEventData> events)
    {
        double eventPotentialTotal = 0.0;
        // Probability of picking each event and their corresponding signature
        foreach (var eventData in events)
        {
            if (eventData.CNEventPars.Size > 0 && McParams.IncludeSize)
            {
                eventPotentialTotal += Math.Log(eventData.GetProb());
            }
            eventPotentialTotal += Math.Log(eventData.CNEventPars.Prob);
        }
        return eventPotentialTotal;
    }

    public double CalculatePotential(double proposedFitness, List<BaseEventData> events)
    {
        double fitnessPotential = McParams.ThetaFitness * proposedFitness;
        return McParams.IncludeProb ? GetEventPotential(events) + fitnessPotential : fitnessPotential;
    }

    public double GetFitness(Karyotype kar, List<BaseEventData> events)
    {
        // Probability of picking each event and their corresponding signature
        foreach (var eventData in events)
            eventData.ApplyEvent(kar);
        return kar.UpdateFitness(GenRef, FitnessParams);
    }

    private BaseEventData GetNewEvent(Sample sample, Karyotype kar)
    {
        var cnEventP = Rnd.PickRndElem(sample.EventPars);
        var eventData = Sampling.GenerateCNEventData(Rnd, kar, cnEventP);
        return eventData ?? throw new Exception("Failed to generate new event data.");
    }

    private List<BaseEventData> Evolve(Sample sample, Karyotype kar)
    {
        var currentEvents = new List<BaseEventData>();
        var currentFitness = Fitness.Calculate(new Karyotype(kar), GenRef, FitnessParams);
        var currentPotential = CalculatePotential(currentFitness, currentEvents);

        var bestFitness = currentFitness;
        var bestEvents = new List<BaseEventData>(currentEvents);

        var fitDict = new Dictionary<int, (int, double)>{};

        fitDict[-1] = (0, currentFitness);

        for (int i = 0; i < McParams.NumSamplesTotal; i++)
        {
            // Generate a new event and correspondingly add to list
            var newEvent = GetNewEvent(sample, new Karyotype(kar));
            var proposedEvents = new List<BaseEventData>(currentEvents) { newEvent };
            var proposedFitness = GetFitness(new Karyotype(kar), proposedEvents);
            var proposedPotential = CalculatePotential(proposedFitness, proposedEvents);
            var acceptProb = proposedPotential - currentPotential;
            if (acceptProb >= Math.Log(Rnd.NextDouble()))
            {
                currentPotential = proposedPotential;
                currentEvents = proposedEvents;
                if (proposedFitness > bestFitness)
                {
                    bestFitness = proposedFitness;
                    bestEvents = proposedEvents;
                }
                fitDict[i] = (currentEvents.Count, proposedFitness);
                // Apply the new event to the clone
                newEvent.ApplyEvent(kar);
            }
        }
        if (McParams.PrintFitnesses)
        {
            FileIO.WriteFitnesses(fitDict);
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

            var bestEvents = Evolve(sample, childKar);

            for (int mutNo = 0; mutNo < bestEvents.Count; mutNo++)
            {
                Console.Write($"\rSample {sample.SampleId}. Clone {Counter}/{clones.Count}. Event {mutNo + 1}/{bestEvents.Count}.");
                var eventData = bestEvents[mutNo];
                eventData.ApplyEvent(childKar);
                double newFitness = childKar.UpdateFitness(GenRef, FitnessParams);
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
