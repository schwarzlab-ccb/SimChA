using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;
using EDists = Extreme.Statistics.Distributions;

namespace SimChA.Simulation;

// TODO: Split into two classes (one for random sampling, one for applying MCMC).
public class Simulator
{
    private readonly Random _rnd;
    private readonly FitnessParams _fitness;
    private readonly Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> _geneLists;
    private int _counter;

    public Simulator(
        Random rnd,
        FitnessParams fitnessParams,
        Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> geneLists)
    {
        _rnd = rnd;
        _fitness = fitnessParams;
        _geneLists = geneLists;
    }
    
    public void SampleEvents(Sample sample)
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
    
    public void SampleEvents(Sample sample, MCParams mcParams)
    {
        if (sample.EventPars == null || !sample.EventPars.Any())
        {
            throw new Exception("No events to sample from.");
        }
        _counter = 1;
        var (root, childLoopUp) = CloneComp.CreateLookUp(sample.Clones);
        sample.Kars[root.CloneId] = new Karyotype(sample.SexXX);
        MCSampleCNEventsRec(sample, root, childLoopUp, mcParams, 1);
    }

    private void ApplyCNEventsRec(Sample sample, CloneIn node, Dictionary<int, List<CloneIn>> clones, int eventCount)
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
                var eventP = Sampling.PickRndElem(_rnd, sample.EventPars);
                var eventData = Sampling.GenerateCNEventData(_rnd, childKar, eventP);
                eventData.ApplyEvent(childKar);
                double newFitness = childKar.UpdateFitness(_geneLists, _fitness);
                double dFit = newFitness - oldFitness;
                var abberation = new CNEventDesc(eventP.Type, eventCount + mutNo,eventData.ToString(), dFit, newFitness);
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
            select new Sample(profile.Key, profile.Value.SexXX, clones, new List<CNEventPars>())).ToList();

    public List<BaseEventData> InitEvents(Karyotype kar, int nMutations, List<CNEventPars> cnEventPs)
    {
        var eventPs = Enumerable.Range(0, nMutations).Select(_ => Sampling.PickRndElem(_rnd, cnEventPs));
        return eventPs.Select(e => Sampling.GenerateCNEventData(_rnd, kar, e)).ToList();
    }

    // The conditional (un-normalized) probability of this set of events occuring, 
    // given the individual events and the signature
    public double Potential(MCParams mcParams, Karyotype kar, double targetFit, List<BaseEventData> events, ref bool thresholdAccept)
    {
        double eventPotentialTotal = 1.0;

        // Probability of picking each event and their corresponding signature
        // (ignore normalization, divides out when we do the accept/reject)
        foreach (var eventData in events)
        {
            eventData.ApplyEvent(kar);
            eventPotentialTotal *= eventData.CNEventPars.Prob;
        }

        double newFitness = kar.UpdateFitness(_geneLists, _fitness);
        double dFit = newFitness - targetFit;
        thresholdAccept = Math.Abs(dFit / targetFit) < mcParams.ThresholdFit;

        // Fitness potential is an exponential - exp[-theta * |fit - mean_fit|]
        double fitnessPotential = Math.Exp(-mcParams.ThetaFitness * Math.Abs(dFit));
        // Gaussian form
        //double fitnessPotential = Math.Exp(-McParams.ThetaFitness*Math.Pow(dFit,2));

        return fitnessPotential * eventPotentialTotal;
    }
    
    private void MCSampleCNEventsRec(Sample sample, CloneIn node, IReadOnlyDictionary<int, List<CloneIn>> clones, MCParams mcParams, int eventCount)
    {
        foreach (var child in clones[node.CloneId])
        {
            var childKar = new Karyotype(sample.Kars[node.CloneId]);
            sample.Kars[child.CloneId] = childKar;
            var childEvs = new List<CNEventDesc>();
            sample.EventDescs[child.CloneId] = childEvs;
            
            // Skip any clones that don't have any mutational distance from their
            // parent
            if (child.Distance <= 0)
                continue;

            // Initialize all the relevant quantities
            double oldFitness = childKar.FitnessVal;

            // Parameters needed for the MH algorithm
            float alterEventStart = 0.5f;
            float alterEventLength = 0.5f;
            // The tracker to accept sets of events if they are within the threshold tolerance from
            // the target fitness
            bool thresholdAccept = false;

            // Generate a starting set of mutations and its potential
            var currentEventProps = InitEvents(childKar, child.Distance, sample.EventPars);
            double currentPotential = Potential(mcParams, new Karyotype(childKar), child.FitnessTarget, currentEventProps, ref thresholdAccept);

            // Now we perform the Metropolis-Hastings algorithm
            // and sample a set of events that give the closest agreement with
            // fitness given by SMITH
            for (int i = 0; i < mcParams.NumSamplesTotal; i++)
            {
                // Reset the automatic acceptance
                thresholdAccept = false;
                var proposedEventProps = currentEventProps.ToList();
                // Select a random CNEventPars to modify
                int index = _rnd.Next(proposedEventProps.Count);
                // Choose whether to swap the event entirely
                if (_rnd.NextDouble() < mcParams.SwapEventP)
                {
                    // Get the new signature and the corresponding event
                    var cnEventP = Sampling.PickRndElem(_rnd, sample.EventPars);
                    proposedEventProps[index] = Sampling.GenerateCNEventData(_rnd, childKar, cnEventP);
                }
                // Otherwise we modify some quantity of the event, but keep the event itself the same
                else
                {
                    // Keep the event type the same, but redo all parameters:
                    var cnEventP = proposedEventProps[index].CNEventPars;
                    proposedEventProps[index] = Sampling.GenerateCNEventData(_rnd, childKar, cnEventP);
                }
                // With the newly selected event, we need to calculate the new
                // fitness of the clone
                double proposalPotential = Potential(mcParams,new Karyotype(childKar), child.FitnessTarget, proposedEventProps, ref thresholdAccept);
                double acceptProb = proposalPotential / currentPotential;
                if (acceptProb >= _rnd.NextDouble())
                {
                    currentPotential = proposalPotential;
                    currentEventProps = proposedEventProps;
                    // Break out of the sampling if we have reached the threshold
                    // and have reached the minimum number of samples required
                    if (thresholdAccept && i > mcParams.NumSamplesMin)
                        break;
                }
            }

            // Finalize the mutated karyotype by applying the best-fit set of events
            // and move onto the next clone
            // Copy its parent
            for (int mutNo = 0; mutNo < currentEventProps.Count; mutNo++)
            {
                Console.Write($"\rClone {_counter}/{clones.Count - 1}. Event {mutNo + 1}/{child.Distance}.");
                var eventData = currentEventProps[mutNo];
                eventData.ApplyEvent(childKar);
                double newFitness = childKar.UpdateFitness(_geneLists, _fitness);
                double dFit = newFitness - oldFitness;
                
                var abberation = new CNEventDesc(eventData.EventType, eventCount + mutNo, eventData.ToString(), dFit, newFitness);
                childEvs.Add(abberation);
            }
            _counter++;            
            if (child.CloneId != node.CloneId)
            {
                MCSampleCNEventsRec(sample, child, clones, mcParams, eventCount + child.Distance);
            }
        }
    }
}