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

    protected static BaseEventData CreatePassEvent() 
        => new(new CNEventPars(CNEventType.Pass, 1));

    public List<Sample> Simulate(CTreeNode root, List<CTreeNode> cloneTree, List<Signature> sigs)
    {
        var (cnEventPs, mixture) = Factory.MixSignatures(sigs);
        var res = new List<Sample>();
        var sex = Sampling.GetSex(Rnd, SimParams.Sex);
        var rootKar = new Karyotype(GenRef, sex);
        if (SimParams.TetraploidStart)
        {
            rootKar.ApplyWGD();
        }
        rootKar.UpdateFitness(GenRef, FitParams);
        ApplyCNEventsRec(root, cloneTree, cnEventPs, mixture, res, rootKar, 0);
        return res;
    }

    protected int SampleEventCount(CTreeNode node)
    {
        double events = node.Distance > 0 ? node.Distance : SimParams.RateMean;
        return Sampling.SampleDistInt(Rnd, SimParams.RateDist, events);
    }
    
    protected double SampleFit(CTreeNode node)
        => node.Fitness > 0 ? node.Fitness : Sampling.SampleDist(Rnd, SimParams.FitDist, SimParams.FitMean);

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
            Console.Write($"\rSample {cnChild.CloneId}. Event {evNo}/{eventCount}.".PadRight(80));
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
            var newClone = new Sample(parent.CloneId, child.CloneId, childKar, childEvs, mixture);
            sampleList.Add(newClone);
            
            if (child.CloneId != parent.CloneId)
            {
                ApplyCNEventsRec(child, cloneTree, cnEventPs, mixture, sampleList, childKar, mutDepth + childEvs.Count);
            }
        }
    }
}
