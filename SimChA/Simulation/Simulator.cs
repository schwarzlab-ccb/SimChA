using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;
using EDists = Extreme.Statistics.Distributions;

namespace SimChA.Simulation;

public class Simulator
{
    public readonly Random _rnd;
    public readonly FitnessParams _fitness;
    public readonly Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> _geneLists;
    public int _counter;

    public Simulator(
        Random rnd,
        FitnessParams fitnessParams,
        Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> geneLists)
    {
        _rnd = rnd;
        _fitness = fitnessParams;
        _geneLists = geneLists;
    }

    public virtual void SampleEvents(Sample sample)
    {
        if (sample.EventPars == null || !sample.EventPars.Any())
        {
            throw new Exception("No events to sample from.");
        }
        _counter = 1;
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
            var childEvs = new List<CNEventDesc>();
            sample.EventDescs[child.CloneId] = childEvs;
            double oldFitness = childKar.FitnessVal;
            for (int mutNo = 0; mutNo < child.Distance; mutNo++)
            {
                Console.Write($"\rClone {_counter}/{clones.Count}. Event {mutNo + 1}/{child.Distance}.");
                var eventP = _rnd.PickRndElem(sample.EventPars);
                var eventData = Sampling.GenerateCNEventData(_rnd, childKar, eventP);
                eventData.ApplyEvent(childKar);
                double newFitness = childKar.UpdateFitness(_geneLists, _fitness);
                double dFit = newFitness - oldFitness;
                var abberation = new CNEventDesc(eventP.Type, eventCount + mutNo, eventData.ToString(), dFit, newFitness);
                childEvs.Add(abberation);
                oldFitness = newFitness;
            }
            _counter++;
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
        var eventPs = Enumerable.Range(0, nMutations).Select(_ => _rnd.PickRndElem(cnEventPs));
        return eventPs.Select(e => Sampling.GenerateCNEventData(_rnd, kar, e)).ToList();
    }
}