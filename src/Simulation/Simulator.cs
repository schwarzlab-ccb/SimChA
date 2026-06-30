using SimChA.Computation;
using SimChA.Data;
using SimChA.EventData;
using SimChA.IO;

namespace SimChA.Simulation;

public class Simulator(Random rnd, RefGen refGen, SimParams simParams, FitParams fitParams)
{
    protected Random Rnd  { get; } = rnd;
    protected RefGen RefGen  { get; } = refGen;
    protected FitParams FitParams  { get; } = fitParams;
    protected SimParams SimParams { get; } = simParams;

    protected static BaseEventData CreateSkipEvent() 
        => new(new CNEventPars(CNEventType.Skip, 1));

    public List<Sample> Simulate(CTreeNode root, List<CTreeNode> cloneTree, List<Signature> sigs)
    {
        var (cnEventPs, sampleMixture) = Factory.MixSignatures(Rnd, sigs, SimParams.Mixture);
        var res = new List<Sample>();
        var sex = SimParams.AutosomesOnly ? SexType.Any : Sampling.GetSex(Rnd, SimParams.Sex);
        var rootKar = new Karyotype(RefGen, sex);
        if (SimParams.TetraploidStart)
        {
            rootKar.ApplyWGD();
        }
        rootKar.UpdateFitness(RefGen, FitParams);
        ApplyCNEventsRec(root, cloneTree, cnEventPs, sampleMixture, res, rootKar, 0);
        return res;
    }

    protected int SampleEventCount(CTreeNode node) 
        => node.Distance >= 0 ? node.Distance : Sampling.SampleDiscDist(Rnd, SimParams.RateDist, SimParams.RateMean);

    protected static (List<string> regionsGained, List<string> regionsLost) GetKaryotypeDiff(
        Dictionary<int, List<string>> regionsBefore, Dictionary<int, List<string>> regionsAfter)
    {
        var regionsGained = new List<string>();
        var regionsLost = new List<string>();
        
        // Regions from whole contigs that were gained
        foreach (var id in regionsAfter.Keys.Where(id => !regionsBefore.ContainsKey(id)))
            regionsGained.AddRange(regionsAfter[id]);
        
        // Regions from whole contigs that were lost
        foreach (var id in regionsBefore.Keys.Where(id => !regionsAfter.ContainsKey(id)))
            regionsLost.AddRange(regionsBefore[id]);
        
        // Regions changed within contigs that exist both before and after
        foreach (var id in regionsBefore.Keys.Intersect(regionsAfter.Keys))
        {
            var before = regionsBefore[id];
            var after = regionsAfter[id];
            var beforeCounts = before.GroupBy(r => r).ToDictionary(g => g.Key, g => g.Count());
            var afterCounts = after.GroupBy(r => r).ToDictionary(g => g.Key, g => g.Count());
            foreach (var key in afterCounts.Keys)
            {
                int diff = afterCounts[key] - (beforeCounts.TryGetValue(key, out int bc) ? bc : 0);
                for (int i = 0; i < diff; i++) regionsGained.Add(key);
            }
            foreach (var key in beforeCounts.Keys)
            {
                int diff = beforeCounts[key] - (afterCounts.TryGetValue(key, out int ac) ? ac : 0);
                for (int i = 0; i < diff; i++) regionsLost.Add(key);
            }
        }
        return (regionsGained, regionsLost);
    }

    protected static (string regionsGained, string regionsLost) CalcKaryotypeDiff(
        Karyotype beforeKar,
        Karyotype afterKar)
    {
        if (!CNEventDesc.PrintDelta)
        {
            return ("", "");
        }

        var regionsBefore = beforeKar.GetRegionDescriptions();
        var regionsAfter = afterKar.GetRegionDescriptions();
        var (regionsGained, regionsLost) = GetKaryotypeDiff(regionsBefore, regionsAfter);
        return (
            "[" + string.Join(",", regionsGained) + "]",
            "[" + string.Join(",", regionsLost) + "]");
    }

    protected virtual (Karyotype childKar, List<CNEventDesc> childEvs) SampleEvents(
        Karyotype parentKar, 
        CTreeNode cnChild, 
        List<CNEventPars> cnEventPs, 
        int mutDepth)
    {
        var currentKar = new Karyotype(parentKar);
        var childEvs = new List<CNEventDesc>();
        int eventCount = SampleEventCount(cnChild);
        for (int evNo = 1; evNo <= eventCount; evNo++)
        {
            Console.Write($"Sample {cnChild.CloneId}. Event {evNo}/{eventCount}.".PadRight(80) + "\r");
            var eventP = Rnd.PickRndElem(cnEventPs);
            var eventData = Sampling.GenerateCNEventData(Rnd, currentKar, eventP) ?? CreateSkipEvent();
            var childKar = new Karyotype(currentKar);
            eventData.ApplyEvent(childKar);
            (string gainedStr, string lostStr) = CalcKaryotypeDiff(parentKar, childKar);
            string karStr = CNEventDesc.PrintKaryotype ? childKar.ToString() : "";
            var newEv = new CNEventDesc(eventData, mutDepth + evNo, Signature: eventP.Signature,
                RegionsGained: gainedStr, RegionsLost: lostStr, Karyotype: karStr);
            childEvs.Add(newEv);
            currentKar = childKar;
        } 
        return (currentKar, childEvs);
    }

    // Maximum number of times a single sample's event selection is restarted when it exceeds
    // SimParams.MaxWGD before the run is aborted to avoid an unbounded loop.
    private const int MaxWgdRestarts = 10000;

    // Wraps SampleEvents with the SimParams.MaxWGD limit: if the generated sample contains more
    // whole-genome doublings than allowed, the entire event selection for the sample is restarted.
    private (Karyotype childKar, List<CNEventDesc> childEvs) SampleEventsLimited(
        Karyotype parentKar,
        CTreeNode child,
        List<CNEventPars> cnEventPs,
        int mutDepth)
    {
        for (int attempt = 1; ; attempt++)
        {
            var (childKar, childEvs) = SampleEvents(parentKar, child, cnEventPs, mutDepth);
            if (SimParams.MaxWGD < 0)
            {
                return (childKar, childEvs);
            }
            int wgdCount = childEvs.Count(e => e.EventData.EventType == CNEventType.WholeGenomeDoubling);
            if (wgdCount <= SimParams.MaxWGD)
            {
                return (childKar, childEvs);
            }
            if (attempt >= MaxWgdRestarts)
            {
                throw new Exception(
                    $"Sample {child.CloneId} exceeded MaxWGD ({SimParams.MaxWGD}) on every one of " +
                    $"{MaxWgdRestarts} restarts. Increase SimParams.MaxWGD or lower the WGD probability.");
            }
            Console.Write(
                $"\rSample {child.CloneId} produced {wgdCount} WGDs (> MaxWGD {SimParams.MaxWGD}), restarting (attempt {attempt}).".PadRight(80));
        }
    }

    private void ApplyCNEventsRec(
        CTreeNode parent,
        List<CTreeNode> cloneTree,
        List<CNEventPars> cnEventPs,
        Dictionary<string, double> mixture,
        List<Sample> sampleList,
        Karyotype parentKar,
        int mutDepth)
    {
        var children = cloneTree.Where(c => c.ParentId == parent.CloneId).ToList();
        foreach (var child in children)
        {
            var (childKar, childEvs) = SampleEventsLimited(parentKar, child, cnEventPs, mutDepth);
            var newClone = new Sample(child.CloneId,parent.CloneId, childKar, childEvs, mixture, child.Fitness);
            sampleList.Add(newClone);
            
            if (child.CloneId != parent.CloneId)
            {
                ApplyCNEventsRec(child, cloneTree, cnEventPs, mixture, sampleList, childKar, mutDepth + childEvs.Count);
            }
        }
    }
}
