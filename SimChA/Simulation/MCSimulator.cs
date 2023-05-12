using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;
using EDists = Extreme.Statistics.Distributions;

namespace SimChA.Simulation;

// TODO: Split into two classes (one for random sampling, one for applying MCMC).
class MCSimulator : Simulator
{
    private MCParams _mcParams;
    public MCSimulator(Random rnd,
        FitnessParams fitnessParams, Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> geneLists,
        MCParams mCParams) : base(rnd, fitnessParams, geneLists)
    {
        _mcParams = mCParams;

    }
    public override void SampleEvents(Sample sample)
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
    public (double potential, bool accept) Potential(Karyotype kar, double targetFit, List<BaseEventData> events)
    {
        double eventPotentialTotal = 0.0;

        // Probability of picking each event and their corresponding signature
        foreach (var eventData in events)
        {
            eventData.ApplyEvent(kar);
            eventPotentialTotal += Math.Log(eventData.CNEventPars.Prob);
        }
        double dFit = kar.UpdateFitness(_geneLists, _fitness) - targetFit;
        // Variable to immediately quit the MC Sampling if we've reached enough accuracy
        bool accept = Math.Abs(dFit / targetFit) < _mcParams.ThresholdFit;
        // Fitness potential is an exponential - exp[-theta * |fit - mean_fit|]
        double fitnessPotential = Math.Exp(-_mcParams.ThetaFitness * Math.Abs(dFit));

        double potential = eventPotentialTotal + Math.Log(fitnessPotential);

        return (potential, accept);
    }
    private List<BaseEventData> GetNewProposal(Sample sample, Karyotype kar, List<BaseEventData> oldEvents)
    {
        var proposedEvents = oldEvents.ToList();
        // Select a random CNEventPars to modify
        int index = _rnd.Next(proposedEvents.Count);
        // Choose whether to swap the event entirely
        if (_rnd.NextDouble() < _mcParams.SwapEventP)
        {
            // Get the new signature and the corresponding event
            var cnEventP = _rnd.PickRndElem(sample.EventPars);
            proposedEvents[index] = Sampling.GenerateCNEventData(_rnd, kar, cnEventP);
        }
        // Otherwise we modify some quantity of the event, but keep the event itself the same
        else
        {
            // Keep the event type the same, but redo the parameters:
            var cnEventP = proposedEvents[index].CNEventPars;
            proposedEvents[index] = Sampling.GenerateCNEventData(_rnd, kar, cnEventP);
        }
        return proposedEvents;
    }
    private List<BaseEventData> GetBestEvents(Sample sample, Karyotype kar, int nEvents, double targetFitness){
        // Generate a starting set of mutations and its potential
        var currentEvents = InitEvents(kar, nEvents, sample.EventPars);
        double currentPotential = Potential(new Karyotype(kar), targetFitness, currentEvents).potential;

        // Now we perform the Metropolis-Hastings algorithm
        // and sample a set of events that give the closest agreement with fitness given by SMITH
        for (int i = 0; i < _mcParams.NumSamplesTotal; i++)
        {
            var proposedEvents = GetNewProposal(sample, kar, currentEvents);
            // Calculate the new fitness of the proposed set of events on the clone
            (double proposalPotential, bool thresholdAccept) = Potential(new Karyotype(kar), targetFitness, proposedEvents);
            double acceptProb = proposalPotential - currentPotential;
            if (acceptProb >= Math.Log(_rnd.NextDouble()))
            {
                currentPotential = proposalPotential;
                currentEvents = proposedEvents;
                // Break out of the sampling if we have reached the threshold
                // and have reached the minimum number of samples required
                if (thresholdAccept && i > _mcParams.NumSamplesMin)
                    break;
            }
        }
        return currentEvents;
    }
    private void ApplyCNEventsRec(Sample sample, CloneIn node, IReadOnlyDictionary<int, 
        List<CloneIn>> clones, int eventCount)
    {
        foreach (var child in clones[node.CloneId])
        {
            var childKar = new Karyotype(sample.Kars[node.CloneId]);
            sample.Kars[child.CloneId] = childKar;
            var childEvs = new List<CNEventDesc>();
            sample.EventDescs[child.CloneId] = childEvs;
            
            if (child.Distance > 0)
            {
                double oldFitness = childKar.FitnessVal;
                // Perform the MC Sampling to generate a set of events
                var bestEvents = GetBestEvents(sample, childKar, child.Distance,child.FitnessTarget);

                // Finalize the mutated karyotype by applying the best-fit set of events
                for (int mutNo = 0; mutNo < bestEvents.Count; mutNo++)
                {
                    Console.Write($"\rClone {_counter}/{clones.Count}. Event {mutNo + 1}/{child.Distance}.");
                    var eventData = bestEvents[mutNo];
                    eventData.ApplyEvent(childKar);
                    double newFitness = childKar.UpdateFitness(_geneLists, _fitness);
                    double dFit = newFitness - oldFitness;
                    var abberation = new CNEventDesc(eventData.EventType, eventCount + mutNo, eventData.ToString(), dFit,
                        newFitness);
                    childEvs.Add(abberation);
                }
                _counter++;
                if (child.CloneId != node.CloneId)
                {
                    ApplyCNEventsRec(sample, child, clones, eventCount + child.Distance);
                }
            }
        }
    }
}