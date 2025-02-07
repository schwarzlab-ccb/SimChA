using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;
using SimChA.IO;

namespace SimChA.Simulation;

public class Simulator
{
    protected Random Rnd  { get; }
    public GenRef GenRef  { get; }
    public FitnessParams FitnessParams  { get; }
    
    protected int Counter;

    public Simulator(Random rnd, GenRef genRef, FitnessParams fitnessParams)
    {
        Rnd = rnd;
        GenRef = genRef;
        FitnessParams = fitnessParams;
    }

    public virtual void SampleEvents(Sample sample)
    {
        if (sample.EventPars == null || !sample.EventPars.Any())
        {
            throw new Exception("No events to sample from.");
        }
        Counter = 1;
        var (root, childLookUp) = CloneComp.CreateLookUp(sample.Clones);
        sample.Kars[root.CloneId] = new Karyotype(GenRef, sample.Sex);
        // TODO: Check if 1 is the correct number of events!
        ApplyCNEventsRec(sample, root, childLookUp, 1);
    }
    
    private void ApplyCNEventsRec(Sample sample, CloneIn node, IReadOnlyDictionary<string, List<CloneIn>> clones, int eventCount)
    {
        foreach (var child in clones[node.CloneId])
        {
            var childKar = new Karyotype(sample.Kars[node.CloneId]);
            sample.Kars[child.CloneId] = childKar;
            var childEvs = new List<CNEventDesc>();
            sample.EventDescs[child.CloneId] = childEvs;
            for (int mutNo = 0; mutNo < child.Distance; mutNo++)
            {
                Console.Write($"\rSample {sample.SampleId}. Clone {Counter}/{clones.Count}. Event {mutNo + 1}/{child.Distance}.".PadRight(80));
                var eventP = Rnd.PickRndElem(sample.EventPars);
                var eventData = Sampling.GenerateCNEventData(Rnd, childKar, eventP);
                // TODO: we should log this somewhere for the user to know that we didn't sample the exact number of events
                if (eventData == null)
                    continue;
                eventData.ApplyEvent(childKar);
                var abberation = new CNEventDesc(eventP.Type, eventCount + mutNo, eventData.ToString());
                childEvs.Add(abberation);
            }
            Counter++;
            if (child.CloneId != node.CloneId)
            {
                ApplyCNEventsRec(sample, child, clones, eventCount + child.Distance);
            }
        }
    }

    public List<BaseEventData> InitEvents(Karyotype kar, int nMutations, List<CNEventPars> cnEventPs)
    {
        var eventPs = Enumerable.Range(0, nMutations).Select(_ => Rnd.PickRndElem(cnEventPs));
        return eventPs.Select(
            e =>
            {
                var newEventD = Sampling.GenerateCNEventData(Rnd, kar, e) 
                                ?? throw new Exception("Failed to initialize event data.");
                return newEventD;
            }
        ).ToList();
    }
    
    public List<(double fitness, int eventCount)> FitnessListFromSamples(SimParams simParams, Dictionary<string, Karyotype> profiles, Dictionary<string, int> eventCounts)
    {
        var output = new List<(double fitness, int eventCount)>();
        var samples = Converters.SamplesFromProfiles(profiles);
        foreach (var sample in samples)
        {
            int total = sample.Clones.Count;
            foreach (var clone in sample.Clones)
            {
                sample.CloneStats[clone.CloneId] = CNProfile.GetCloneStats(sample, clone, GenRef, simParams.Fitness, sample.Kars);
                output.Add((sample.CloneStats[clone.CloneId].Fitness, eventCounts[sample.SampleId]));
            }
        }
        return output;
    }
}
