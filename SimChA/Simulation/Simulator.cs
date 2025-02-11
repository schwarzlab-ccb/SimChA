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

    public List<Sample> Simulate(CTreeNode root, List<CTreeNode> cloneTree, List<Signature> sigs)
    {
        var (cnEventPs, mixture) = Converters.PropagateSigs(sigs);
        var res = new List<Sample>();
        var sex = Sampling.GetSex(Rnd, SimParams.Sex);
        var rootKar = new Karyotype(GenRef, sex);
        if (SimParams.TetraploidStart)
        {
            rootKar.ApplyWGD();
            rootKar.UpdateFitness(GenRef, FitParams);
        }
        ApplyCNEventsRec(root, cloneTree, cnEventPs, mixture, res, rootKar, 1);
        return res;
    }
    
    protected virtual void ApplyCNEventsRec(
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
            var childKar = new Karyotype(parentKar);
            var childEvs = new List<CNEventDesc>();
            int distance = child.Distance > 0
                ? child.Distance
                : Sampling.SampleDistInt(Rnd, SimParams.RateDist, SimParams.RateMean);
            for (int mutNo = 0; mutNo < distance; mutNo++)
            {
                Console.Write($"\rSample {child.CloneId}. Event {mutNo + 1}/{child.Distance}.".PadRight(80));
                var eventP = Rnd.PickRndElem(cnEventPs);
                var eventData = Sampling.GenerateCNEventData(Rnd, childKar, eventP);
                // TODO: we should log this somewhere for the user to know that we didn't sample the exact number of events
                if (eventData == null)
                    continue;
                eventData.ApplyEvent(childKar);
                var abberation = new CNEventDesc(eventP.Type, mutDepth + mutNo, eventData.ToString());
                childEvs.Add(abberation);
            }

            var newClone = new Sample(parent.CloneId, child.CloneId, childKar, childEvs, mixture);
            sampleList.Add(newClone);
            if (child.CloneId != parent.CloneId)
            {
                ApplyCNEventsRec(child, cloneTree, cnEventPs, mixture, sampleList, childKar, mutDepth + childEvs.Count);
            }
        }
    }
}
