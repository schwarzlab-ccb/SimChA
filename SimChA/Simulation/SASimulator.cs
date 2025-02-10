using SimChA.Computation;
using SimChA.Data;
using SimChA.EventData;
using SimChA.IO;
namespace SimChA.Simulation;

public class SASimulator : Simulator
{
    private EvoParams EvoParams { get; }

    public SASimulator(Random rnd, GenRef genRef, SimParams simParams, FitParams fitParams, EvoParams evoParams) 
        : base(rnd, genRef, simParams, fitParams)
    {
        EvoParams = evoParams;
    }

    public override void Simulate(Sample sample)
    {
        if (sample.EventPars == null || sample.EventPars.Count == 0)
        {
            throw new Exception("No events to sample from.");
        }
        var (root, childLookUp) = CloneComp.CreateLookUp(sample.Clones);
        sample.Kars[root.CloneId] = new Karyotype(GenRef, sample.Sex);
        // Start with the tetraploid state
        if (SimParams.TetraploidStart)
        {
            sample.Kars[root.CloneId].ApplyWGD();
            sample.Kars[root.CloneId].UpdateFitness(GenRef, FitParams);
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

    private BaseEventData GetNewEvent(List<CNEventPars> cnEventPars, Karyotype kar, double currentFitness)
    {
        List<(double, BaseEventData)> events = new();
        for (int tryNo = 0; tryNo <= EvoParams.MaxTries; tryNo++)
        {
            var cnEventP = Rnd.PickRndElem(cnEventPars);
            var eventData = Sampling.GenerateCNEventData(Rnd, kar, cnEventP);
            if (eventData == null)
            {
                continue;
            }
            var proposedKar = new Karyotype(kar);
            eventData.ApplyEvent(proposedKar);
            double proposedFitness = proposedKar.UpdateFitness(GenRef, FitParams);
            if (Math.Exp(proposedFitness - currentFitness - EvoParams.Acceptance) >= Rnd.NextDouble())
            {
                return eventData;
            }
            events.Add((proposedFitness, eventData));
        }
        if (events.Count == 0)
        {
            throw new Exception($"Could not generate a new event. In {EvoParams.MaxTries} tries.");
        }

        var best = events.MaxBy(e => e.Item1);
        return best.Item2;
    }

    bool MetBreakCondition(int time, int nEventsLeft)
        => !EvoParams.EvolveInTime ? nEventsLeft == 0 : time >= 1; 

    double GetMutRate(Karyotype kar)
        => EvoParams.EventRate * CNProfile.CalcPloidy(kar, GenRef);

    bool DidMutate(Karyotype kar)
        => !EvoParams.EvolveInTime || Rnd.NextDouble() < 1 - Math.Exp(-GetMutRate(kar));
    
    private void ApplyEvolutionRec(
        Sample sample,
        CTreeNode node,
        IReadOnlyDictionary<string, List<CTreeNode>> clones,
        int eventCount)
    {
        foreach (var child in clones[node.CloneId])
        {
            int mutCount = child.Distance;
            var childEvs = new List<CNEventDesc>();
            double currentFitness = Fitness.Calculate(sample.Kars[node.CloneId], GenRef, FitParams);
            bool hasWGD = false;
            int time = 0;
            while (!MetBreakCondition(time, mutCount - childEvs.Count))
            {
                time++;
                if (!DidMutate(new Karyotype(sample.Kars[node.CloneId])))
                {
                    continue;
                }
                Console.Write($"\rSample: {sample.SampleId}. " +
                              $"Mutation: {childEvs.Count+1}/{mutCount}.");
                
                var cnEventPars = GetEventPars(sample.EventPars, hasWGD);
                var newEvent = GetNewEvent(cnEventPars, sample.Kars[node.CloneId], currentFitness);

                var newKar = new Karyotype(sample.Kars[node.CloneId]);
                newEvent.ApplyEvent(newKar);
                double proposedFitness = newKar.UpdateFitness(GenRef, FitParams);
                double dFit = proposedFitness - currentFitness;
                var eventTime = EvoParams.EvolveInTime ? time : childEvs.Count + 1;
                var abberation = new CNEventDesc(
                    newEvent.EventType,
                    eventCount + childEvs.Count + 1,
                    newEvent.ToString(),
                    dFit,
                    proposedFitness,
                    eventTime);
                hasWGD |= newEvent.EventType == CNEventType.WholeGenomeDoubling;
                childEvs.Add(abberation);
                sample.Kars[node.CloneId] = newKar;
                currentFitness = proposedFitness;
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