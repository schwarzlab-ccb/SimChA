using SimChA.Computation;
using SimChA.Data;
using SimChA.EventData;
using SimChA.IO;

namespace SimChA.Simulation;

public class Simulator
{
    protected Random Rnd  { get; }
    protected GenRef GenRef  { get; }
    protected FitParams FitParams  { get; }
    protected SampleParams SampleParams { get; }
    
    protected int Counter;

    public Simulator(Random rnd, GenRef genRef, SampleParams sampleParams, FitParams fitParams)
    {
        Rnd = rnd;
        GenRef = genRef;
        FitParams = fitParams;
        SampleParams = sampleParams;
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
            e => Sampling.GenerateCNEventData(Rnd, kar, e) ?? throw new Exception($"Failed to generate event data for {e}.")
        ).ToList();
    }
}
