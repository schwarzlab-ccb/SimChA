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
    protected SimParams SimParams { get; }
    
    public Simulator(Random rnd, GenRef genRef, SimParams simParams, FitParams fitParams)
    {
        Rnd = rnd;
        GenRef = genRef;
        FitParams = fitParams;
        SimParams = simParams;
    }

    public virtual List<Sample> Simulate(Sample sample)
    {
        if (sample.EventPars == null || !sample.EventPars.Any())
        {
            throw new Exception("No events to sample from.");
        }
        var (root, childLookUp) = CloneComp.CreateLookUp(sample.Clones);
        var res = new List<Sample>();
        var rootKar =  new Karyotype(GenRef, sample.Sex);
        res.Add(new Sample(root.CloneId, rootKar, new List<CNEventDesc>()));
        ApplyCNEventsRec(sample, root, childLookUp, res, rootKar, 1);
        return res;
    }
    
    private void ApplyCNEventsRec(
        Sample sample, 
        CTreeNode node, 
        IReadOnlyDictionary<string, List<CTreeNode>> cloneLookUp, 
        List<Sample> clones,
        Karyotype parentKar,
        int eventCount)
    {
        foreach (var child in cloneLookUp[node.CloneId])
        {
            var childKar = new Karyotype(parentKar);
            var childEvs = new List<CNEventDesc>();
            for (int mutNo = 0; mutNo < child.Distance; mutNo++)
            {
                Console.Write($"\rSample {sample.SampleId}. Clone {clones.Count}/{cloneLookUp.Count}. Event {mutNo + 1}/{child.Distance}.".PadRight(80));
                var eventP = Rnd.PickRndElem(sample.EventPars);
                var eventData = Sampling.GenerateCNEventData(Rnd, childKar, eventP);
                // TODO: we should log this somewhere for the user to know that we didn't sample the exact number of events
                if (eventData == null)
                    continue;
                eventData.ApplyEvent(childKar);
                var abberation = new CNEventDesc(eventP.Type, eventCount + mutNo, eventData.ToString());
                childEvs.Add(abberation);
            }

            var newClone = new Sample(child.CloneId, childKar, childEvs);
            clones.Add(newClone);
            if (child.CloneId != node.CloneId)
            {
                ApplyCNEventsRec(sample, child, cloneLookUp, clones, childKar, eventCount + childEvs.Count);
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
