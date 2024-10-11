using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;
using SimChA.IO;
using EDists = Extreme.Statistics.Distributions;

namespace SimChA.Simulation;

public class Evolver
{
    protected readonly FitnessParams FitnessParams;
    protected readonly MCParams McParams;
    protected readonly FileIO FileIO;
    protected readonly Random Rnd;
    protected readonly GenRef GenRef;
    protected int Counter;
    protected List<CNEventPars> EventPars;

    public Evolver(
        Random rnd,
        GenRef genRef,
        FitnessParams fitnessParams, 
        MCParams mCParams,
        FileIO fileIO)
    {
        Rnd = rnd;
        GenRef = genRef;
        FitnessParams = fitnessParams;
        McParams = mCParams;
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
        double fitnessPotential = McParams.ThetaFitness * proposedFitness;
        return fitnessPotential;
    }

    private double InverseEventProb(BaseEventData eventData)
    {
        var inverseEvent = new CNEventType();
        switch (eventData.CNEventPars.Type) {
            case CNEventType.ChromDeletion:
                inverseEvent = CNEventType.ChromDuplication;
                break;
            case CNEventType.ChromDuplication:
                inverseEvent = CNEventType.ChromDeletion;
                break;
            case CNEventType.InternalDeletion:
                inverseEvent = CNEventType.InternalDuplication;
                break;
            case CNEventType.InternalDuplication:
                inverseEvent = CNEventType.InternalDeletion;
                break;
            case CNEventType.TailDeletion:
                inverseEvent = CNEventType.TailDuplication;
                break;
            case CNEventType.TailDuplication:
                inverseEvent = CNEventType.TailDeletion;
                break;
            case CNEventType.ArmDeletion:
                inverseEvent = CNEventType.ArmDuplication;
                break;
            case CNEventType.ArmDuplication:
                inverseEvent = CNEventType.ArmDeletion;
                break;
            case CNEventType.SNV:
            case CNEventType.BreakageFusionBridge:
            case CNEventType.TICycle:
            case CNEventType.TIBridge:
            case CNEventType.TIChain:
            case CNEventType.Rigma:
            case CNEventType.Pyrgo:
            case CNEventType.Chromothripsis:
            case CNEventType.InternalInversion:
            case CNEventType.InvertedDuplication:
            case CNEventType.Chromoplexy:
            case CNEventType.Translocation:
            case CNEventType.WholeGenomeDoubling:
                throw new Exception($"{eventData.CNEventPars.Type} does not currently have an inverse event implemented.");
        }
        var inverseEventData = EventPars.Find(e => e.Type == inverseEvent) ?? throw new Exception($"Could not find inverse event data in the input configuration file for {inverseEvent}.");
        return inverseEventData.Prob;
    }

    private double CalculateTransition(BaseEventData proposedEvent)
    {
        double forwardsProb = Math.Log(proposedEvent.CNEventPars.Prob);
        double backwardsProb = Math.Log(InverseEventProb(proposedEvent));
        return backwardsProb - forwardsProb;
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

    private List<BaseEventData> Evolve(Sample sample, Karyotype kar)
    {
        var currentEvents = new List<BaseEventData>();
        var currentFitness = Fitness.Calculate(new Karyotype(kar), GenRef, FitnessParams);
        var currentPotential = CalculatePotential(currentFitness);

        var bestFitness = currentFitness;

        var fitDict = new Dictionary<int, (int, double)>{};

        fitDict[-1] = (0, currentFitness);

        for (int i = 0; i < McParams.NumSamplesTotal; i++)
        {
            Console.Write($"\rSample {sample.SampleId}. Event {i+1}/{McParams.NumSamplesTotal}.".PadRight(80));
            // Generate a new event and correspondingly add to list
            var newEvent = GetNewEvent(sample, new Karyotype(kar));
            var proposedFitness = GetFitness(new Karyotype(kar), newEvent);
            var proposedPotential = CalculatePotential(proposedFitness);
            var acceptProb = Math.Min(0, proposedPotential - currentPotential + CalculateTransition(newEvent));
            if (acceptProb >= Math.Log(Rnd.NextDouble()))
            {
                currentPotential = proposedPotential;
                currentEvents.Add(newEvent);
                if (proposedFitness > bestFitness)
                {
                    bestFitness = proposedFitness;
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

            /*for (int mutNo = 0; mutNo < bestEvents.Count; mutNo++)
            {
                Console.Write($"\rSample {sample.SampleId}. Clone {Counter}/{clones.Count}. Event {mutNo + 1}/{bestEvents.Count}.");
                var eventData = bestEvents[mutNo];
                eventData.ApplyEvent(childKar);
                double newFitness = childKar.UpdateFitness(GenRef, FitnessParams);
                double dFit = newFitness - oldFitness;
                var abberation = new CNEventDesc(eventData.EventType, eventCount + mutNo, eventData.ToString(), dFit, newFitness);
                childEvs.Add(abberation);
                oldFitness = newFitness;
            }*/

            Counter++;
            if (child.CloneId != node.CloneId)
            {
                ApplyEvolutionRec(sample, child, clones, eventCount + child.Distance);
            }
        }
    }
}
