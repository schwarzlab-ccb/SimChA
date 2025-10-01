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

    protected static BaseEventData CreatePassEvent() 
        => new(new CNEventPars(CNEventType.Pass, 1));

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
    {
        if (node.Distance == 0)
        {
            return 0;
        }
        double events = node.Distance > 0 ? node.Distance : SimParams.RateMean;
        return Sampling.SampleDiscDist(Rnd, SimParams.RateDist, events);
    }
    
    protected virtual (Karyotype childKar, List<CNEventDesc> childEvs) SampleEvents(
        Karyotype parentKar, 
        CTreeNode cnChild, 
        List<CNEventPars> cnEventPs, 
        int mutDepth)
    {
        var childKar = new Karyotype(parentKar);
        var childEvs = new List<CNEventDesc>();
        int eventCount = SampleEventCount(cnChild);
        for (int evNo = 1; evNo <= eventCount; evNo++)
        {
            Console.Write($"Sample {cnChild.CloneId}. Event {evNo}/{eventCount}.".PadRight(80) + "\r");
            var eventP = Rnd.PickRndElem(cnEventPs);
            var eventData = Sampling.GenerateCNEventData(Rnd, childKar, eventP) ?? CreatePassEvent();
            eventData.ApplyEvent(childKar);
            var newEv = new CNEventDesc(eventData, mutDepth + evNo);
            childEvs.Add(newEv);
        } 
        return (childKar, childEvs);
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
            var (childKar, childEvs) = SampleEvents(parentKar, child, cnEventPs, mutDepth);
            var newClone = new Sample(child.CloneId,parent.CloneId, childKar, childEvs, mixture);
            sampleList.Add(newClone);
            
            if (child.CloneId != parent.CloneId)
            {
                ApplyCNEventsRec(child, cloneTree, cnEventPs, mixture, sampleList, childKar, mutDepth + childEvs.Count);
            }
        }
    }
}
