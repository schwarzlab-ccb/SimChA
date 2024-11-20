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
        if (sample.EventPars == null || sample.EventPars.Count == 0)
        {
            throw new Exception("No events to sample from.");
        }
        Counter = 1;
        EventPars = sample.EventPars;
        var (root, childLookUp) = CloneComp.CreateLookUp(sample.Clones);
        sample.Kars[root.CloneId] = new Karyotype(GenRef, sample.Sex);
        ApplyEvolutionRec(sample, root, childLookUp, 1);
    }

    private bool DidMutate(double mean)
        => Rnd.NextDouble() < 1 - Math.Exp(-mean);

    private double CalculateLogAcceptance(double newFitness, double oldFitness, double temperature)
    {
        if (!EvoParams.WithFitness)
        {
            return 0;
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
        return Math.Min(0, fitPart);
    }
    public double GetFitness(Karyotype kar, List<BaseEventData> eventData)
    {
        foreach (var ev in eventData)
        {
            ev.ApplyEvent(kar);
        }
        return kar.UpdateFitness(GenRef, FitnessParams);
    }

    private int GetEventCount(Karyotype kar)
    {
        int nEvents;
        var mu = EvoParams.DynamicMutRate
                ? EvoParams.MutationRate * CNProfile.CalcPloidy(kar, GenRef) / 2.0
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

    private List<CNEventPars> GetModifiedEventPars(List<CNEventPars> pars, Karyotype kar)
    {
        var newPars = new List<CNEventPars>(pars);
        var factor = CNProfile.CalcPloidy(kar, GenRef)/2.0;
        var totalWeight = 0.0;
        foreach (var e in pars)
        {
            var newProb = e.Type switch
            {
                CNEventType.ChromDeletion 
                or CNEventType.ArmDeletion 
                    => Math.Max(0, e.Prob * factor),
                _ => e.Prob,
            };
            totalWeight += newProb;
            newPars[newPars.IndexOf(e)] = e with { Prob = newProb };
        }
        var w = 0.0;
        for (int i = 0; i < newPars.Count; i++)
        {
            newPars[i] = newPars[i] with { Prob = newPars[i].Prob / totalWeight };
            w += newPars[i].Prob;
        }
        return newPars;
    }

    private List<BaseEventData> GetNewEvents(Sample sample, Karyotype kar, int nEvents)
    {
        var sampledEvents = new List<BaseEventData>();
        int iTries = 0;
        var pars = EvoParams.EventCost
                ? GetModifiedEventPars(sample.EventPars, kar)
                : sample.EventPars;
        for (int i = 0; i < nEvents && iTries < EvoParams.MaxTries; )
        {            
            var eventData = GetNewEvent(sample, kar, pars);
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
        /*if (iTries >= EvoParams.MaxTries)
        {
            throw new Exception("Could not generate new events.");
        }*/
        return sampledEvents;
    }

    private BaseEventData? GetNewEvent(Sample sample, Karyotype kar, List<CNEventPars> pars)
    {
        var cnEventP = Rnd.PickRndElem(pars);
        return Sampling.GenerateCNEventData(Rnd, kar, cnEventP);
    }

    private List<BaseEventData> EvolveInContinuousTime(Sample sample, Karyotype kar)
    {
        var currentEvents = new List<BaseEventData>();
        var currentFitness = Fitness.Calculate(new Karyotype(kar), GenRef, FitnessParams);
        var tNow = new List<double>{0.0};
        while (tNow.Last() < EvoParams.MaxTime)
        {
            // Sample the new time for the event
            var u = Rnd.NextDouble();
            var tNew = tNow.Last() - Math.Log(u) / EvoParams.MutationRate;
            if (tNew > EvoParams.MaxTime)
            {
                break;
            }
            tNow.Add(tNew);
            // Generate a new event and correspondingly add to list
            var newEvents = GetNewEvents(sample, new Karyotype(kar), 1);
            if (newEvents.Count != 1)
            {
                continue;
            }
            var ev = newEvents[0];
            var proposedFitness = GetFitness(new Karyotype(kar), newEvents);
            var acceptProb = CalculateLogAcceptance(proposedFitness, currentFitness, EvoParams.Temperature);
            if (acceptProb >= Math.Log(Rnd.NextDouble()))
            {
                currentFitness = proposedFitness;
                currentEvents.Add(ev);
                ev.ApplyEvent(kar);
                kar.UpdateFitness(GenRef, FitnessParams);
            }
        }

        return currentEvents;
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
            // Want to sample a number of events.
            int nEvents = GetEventCount(kar);
            var newEvents = GetNewEvents(sample, new Karyotype(kar), nEvents);
            if (newEvents.Count == 0)
            {
                continue;
            }
            var proposedFitness = GetFitness(new Karyotype(kar), newEvents);
            var acceptProb = CalculateLogAcceptance(proposedFitness, currentFitness, currentTemp);
            if (acceptProb >= Math.Log(Rnd.NextDouble()))
            {
                currentFitness = proposedFitness;
                foreach (var ev in newEvents)
                {
                    currentEvents.Add(ev);
                    ev.ApplyEvent(kar);
                }
                kar.UpdateFitness(GenRef, FitnessParams);
            }
            currentTemp *= EvoParams.SimulatedAnnealing ? EvoParams.CoolingRate : 1.0;
        }
        return currentEvents;
    }

    private int GetNumSteps(int baseNum, Karyotype kar)
    {
        if (EvoParams.MutationRate < 0)
        {
            throw new Exception("Mutation rate must be positive.");
        }
        return (int)Math.Round(baseNum/EvoParams.MutationRate);
    }

    private List<BaseEventData> EvolveInEvents(Sample sample, Karyotype kar, int mutCount)
    {
        var currentEvents = new List<BaseEventData>();
        var currentFitness = Fitness.Calculate(new Karyotype(kar), GenRef, FitnessParams);
        var currentTemp = EvoParams.Temperature;
        
        var nSteps = mutCount;//GetNumSteps(mutCount, kar);
        int i = 0;
        for (; i < nSteps; i++)
        {
            Console.Write($"\rSample {sample.SampleId}. Iteration {i+1}/{nSteps};".PadRight(80));
            // Generate a new event and correspondingly add to list
            int nEvents = GetEventCount(kar);
            var newEvents = GetNewEvents(sample, new Karyotype(kar), nEvents);
            if (newEvents.Count == 0)
            {
                continue;
            }
            var proposedFitness = GetFitness(new Karyotype(kar), newEvents);
            var acceptProb = CalculateLogAcceptance(proposedFitness, currentFitness, currentTemp);
            if (acceptProb >= Math.Log(Rnd.NextDouble()))
            {
                currentFitness = proposedFitness;
                foreach (var ev in newEvents)
                {
                    currentEvents.Add(ev);
                    ev.ApplyEvent(kar);
                }
                kar.UpdateFitness(GenRef, FitnessParams);
            }
            currentTemp *= EvoParams.SimulatedAnnealing ? EvoParams.CoolingRate : 1.0;
        }

        return currentEvents;
    }

    private void ApplyEvolutionRec(Sample sample, CloneIn node, 
        IReadOnlyDictionary<string, List<CloneIn>> clones, 
        int eventCount)
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

            Counter++;
            if (child.CloneId != node.CloneId)
            {
                ApplyEvolutionRec(sample, child, clones, eventCount + child.Distance);
            }
        }
    }
}
