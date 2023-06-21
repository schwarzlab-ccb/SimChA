using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;
using EDists = Extreme.Statistics.Distributions;

namespace SimChA.Simulation;

public class Simulator
{
    protected readonly Random Rnd;
    protected readonly Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> GeneLists;
    protected int Counter;
    protected readonly List<GenContents>? GenContents;

    public Simulator(
        Random rnd,
        Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> geneLists,
        List<GenContents>? genContents = null)
    {
        Rnd = rnd;
        GeneLists = geneLists;
        GenContents = genContents;
    }

    public virtual void SampleEvents(Sample sample)
    {
        if (sample.EventPars == null || !sample.EventPars.Any())
        {
            throw new Exception("No events to sample from.");
        }
        Counter = 1;
        var (root, childLoopUp) = CloneComp.CreateLookUp(sample.Clones);
        sample.Kars[root.CloneId] = new Karyotype(sample.SexXX);
        ApplyCNEventsRec(sample, root, childLoopUp, 1);
    }
    
    private void ApplyCNEventsRec(Sample sample, CloneIn node, IReadOnlyDictionary<int, List<CloneIn>> clones, int eventCount)
    {
        foreach (var child in clones[node.CloneId])
        {
            var childKar = new Karyotype(sample.Kars[node.CloneId]);
            sample.Kars[child.CloneId] = childKar;
            childKar.SetGenContents(GenContents);
            var childEvs = new List<CNEventDesc>();
            sample.EventDescs[child.CloneId] = childEvs;
            for (int mutNo = 0; mutNo < child.Distance; mutNo++)
            {
                Console.Write($"\rSample {sample.SampleId}. Clone {Counter}/{clones.Count}. Event {mutNo + 1}/{child.Distance}.".PadRight(80));
                var eventP = Rnd.PickRndElem(sample.EventPars);
                var eventData = Sampling.GenerateCNEventData(Rnd, childKar, eventP);
                if (eventData == null)
                    return;
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

    public static List<Sample> SamplesFromProfiles(Dictionary<string, Karyotype> profiles)
        => (from profile in profiles
            let clones = new List<CloneIn> { new(0, -1, 0, 0) }
            select new Sample(profile.Key, profile.Value.SexXX, clones, new List<CNEventPars>())
            { Kars = { [0] = profile.Value } }).ToList();

    public List<BaseEventData> InitEvents(Karyotype kar, int nMutations, List<CNEventPars> cnEventPs)
    {
        var eventPs = Enumerable.Range(0, nMutations).Select(_ => Rnd.PickRndElem(cnEventPs));
        return eventPs.Select(
            e =>
            {
                var newEventD = Sampling.GenerateCNEventData(Rnd, kar, e);
                if (newEventD == null)
                {
                    throw new Exception("Failed to initialize event data.");
                }
                return newEventD;
            }
        ).ToList();
    }
}