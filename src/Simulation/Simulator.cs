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
        => node.Distance >= 0 ? node.Distance : Sampling.SampleDiscDist(Rnd, SimParams.RateDist, SimParams.RateMean, SimParams.RateShape);

    protected static (string regionsGained, string regionsLost) CalcKaryotypeDiff(
        Karyotype beforeKar,
        Karyotype afterKar)
    {
        if (!CNEventDesc.PrintDelta)
        {
            return ("", "");
        }

        var gained = new List<string>();
        var lost = new List<string>();
        foreach (var delta in CopyNumbers.DiffKaryotypes(beforeKar, afterKar))
        {
            AppendCopies(gained, lost, "H1", delta.Chrom, delta.Start, delta.End, delta.CNH1);
            AppendCopies(gained, lost, "H2", delta.Chrom, delta.Start, delta.End, delta.CNH2);
        }
        return (
            "[" + string.Join(",", gained) + "]",
            "[" + string.Join(",", lost) + "]");
    }

    // Emits one entry per copy changed on a haplotype: |delta| copies of the interval go to
    // `gained` when the copy number rose and to `lost` when it fell, so the entry count equals
    // the number of copies gained or lost there.
    private static void AppendCopies(
        List<string> gained, List<string> lost,
        string hap, string chrom, long start, long end, int delta)
    {
        if (delta == 0)
        {
            return;
        }
        string desc = $"{hap}:{chrom}[{start}:{end})";
        var target = delta > 0 ? gained : lost;
        for (int i = 0; i < Math.Abs(delta); i++)
        {
            target.Add(desc);
        }
    }

    protected virtual (Karyotype childKar, List<CNEventDesc> childEvs) SampleEvents(
        Karyotype parentKar,
        CTreeNode cnChild,
        List<CNEventPars> cnEventPs,
        int mutDepth,
        int eventCount)
    {
        var currentKar = new Karyotype(parentKar);
        var childEvs = new List<CNEventDesc>();
        for (int evNo = 1; evNo <= eventCount; evNo++)
        {
            Console.Write($"Sample {cnChild.CloneId}. Event {evNo}/{eventCount}.".PadRight(80) + "\r");
            var eventP = Rnd.PickRndElem(cnEventPs);
            var eventData = Sampling.GenerateCNEventData(Rnd, currentKar, eventP) ?? CreateSkipEvent();
            var childKar = new Karyotype(currentKar);
            eventData.ApplyEvent(childKar);
            (string gainedStr, string lostStr) = CalcKaryotypeDiff(currentKar, childKar);
            string karStr = CNEventDesc.PrintKaryotype ? childKar.ToString() : "";
            var newEv = new CNEventDesc(eventData, mutDepth + evNo, Signature: eventP.Signature,
                RegionsGained: gainedStr, RegionsLost: lostStr, Karyotype: karStr);
            childEvs.Add(newEv);
            currentKar = childKar;
        } 
        return (currentKar, childEvs);
    }

    // Absolute cap on the total number of re-simulations for a single sample under the
    // SimParams.MaxWGD limit, guarding against configurations that can never satisfy the cap.
    private const int MaxWgdTotalAttempts = 10000;

    // Wraps SampleEvents with the SimParams.MaxWGD limit. The target event count is drawn once;
    // if the generated sample contains more whole-genome doublings than allowed it is re-simulated
    // with the SAME event count, so the cap does not bias the event-count distribution. The count is
    // only redrawn (and the per-count try budget reset) after SimParams.MaxWgdTries consecutive
    // failures. A fixed tree distance (cnChild.Distance >= 0) cannot be redrawn, so it aborts instead.
    private (Karyotype childKar, List<CNEventDesc> childEvs) SampleEventsLimited(
        Karyotype parentKar,
        CTreeNode child,
        List<CNEventPars> cnEventPs,
        int mutDepth)
    {
        int maxTries = Math.Max(1, SimParams.MaxWgdTries);
        int eventCount = SampleEventCount(child);
        int triesAtCount = 0;
        for (int attempt = 1; ; attempt++)
        {
            var (childKar, childEvs) = SampleEvents(parentKar, child, cnEventPs, mutDepth, eventCount);
            if (SimParams.MaxWGD < 0)
            {
                return (childKar, childEvs);
            }
            int wgdCount = childEvs.Count(e => e.EventData.EventType == CNEventType.WholeGenomeDoubling);
            if (wgdCount <= SimParams.MaxWGD)
            {
                return (childKar, childEvs);
            }
            if (attempt >= MaxWgdTotalAttempts)
            {
                throw new Exception(
                    $"Sample {child.CloneId} exceeded MaxWGD ({SimParams.MaxWGD}) on every one of " +
                    $"{MaxWgdTotalAttempts} re-simulations. Increase SimParams.MaxWGD or lower the WGD probability.");
            }
            if (++triesAtCount >= maxTries)
            {
                if (child.Distance >= 0)
                {
                    throw new Exception(
                        $"Sample {child.CloneId} (fixed event count {eventCount}) exceeded MaxWGD " +
                        $"({SimParams.MaxWGD}) on all {maxTries} re-simulations. A fixed tree distance " +
                        "cannot be redrawn; increase SimParams.MaxWGD or lower the WGD probability.");
                }
                eventCount = SampleEventCount(child);
                triesAtCount = 0;
            }
            Console.Write(
                $"\rSample {child.CloneId} produced {wgdCount} WGDs (> MaxWGD {SimParams.MaxWGD}), re-simulating at count {eventCount} (attempt {attempt}).".PadRight(80));
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
