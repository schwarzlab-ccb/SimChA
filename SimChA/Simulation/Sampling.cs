using Extreme.Statistics.Distributions;
using SimChA.DataTypes;
using SimChA.EventData;
using SimChA.Misc;

namespace SimChA.Simulation;

public static class Sampling
{
    public static long GetNormSeg(Random rnd, long contigLen, double meanFrac) 
        => Math.Max(1, Math.Min((long) Math.Round(contigLen * NormalDistribution.Sample(rnd, meanFrac, meanFrac / 3)), contigLen - 2));
    
    public static long GetExpSeg(Random rnd, long contigLen, long meanLen) 
        => Math.Max(1, Math.Min((long) Math.Round(ExponentialDistribution.Sample(rnd, meanLen)), contigLen - 2));
    
    public static long GetExpSeg(Random rnd, long contigLen, double meanFrac) 
        => Math.Max(1, Math.Min((long) Math.Round(contigLen * ExponentialDistribution.Sample(rnd, meanFrac)), contigLen - 2));

    // Get two positions within the contig (boundaries are excluded)
    public static long GetInternalPos(Random rnd, long contigLen)
        => rnd.NextInt64(1, Math.Max(1, contigLen - 1));
    
    public static int GetFragCount(Random rnd, double mean)
        => (int) Math.Max(1, Math.Round(NormalDistribution.Sample(rnd, mean)));

    // https://ashpublications.org/blood/article/134/Supplement_1/3767/424006/Chromoplexy-and-Chromothripsis-Are-Important
    public static int GetChromoplexySiteCount(Random rnd)
        => rnd.NextSingle() switch
        {
            var n when n < .46 => 3,
            var n when n < .64 => 4,
            var n when n < .74 => 5,
            var n when n < .79 => 6,
            _ => 2
        };

    public static List<long> GetStopsForShards(Random rnd, long contigLen, int shardCount) 
        => Enumerable.Range(0, shardCount - 1)
            .Select(_ => GetInternalPos(rnd, contigLen))
            .Distinct()
            .OrderBy(i => i)
            .ToList();
    
    public static List<double> CreateRandomMixture(Random rnd, double[] concentrations)
        => concentrations.Any() ? new DirichletDistribution(concentrations).Sample(rnd).ToList() : new List<double>();

    public static CNEventP PickRandomEventP(Random rnd, List<CNEventP> eventPs)
    {
        var localEventPs = eventPs.Select(ev => Math.Max(0, ev.Prob) > 0).ToList();
        var probs = localEventPs.Select(ev => ev.Prob).ToList();
        probs = probs.Select(p => p / probs.Sum()).ToList();
        int index = PickRandomIndex(rnd, probs);
        return localEventPs[index];
    }
    
    public static int PickRandomIndex(Random rnd, List<double> probs)
    {
        double val = rnd.NextDouble();
        for (var i = 0; i < probs.Count; i++)
        {
            if (val < probs[i])
            {
                return i;
            }
            val -= probs[i];
        }
        return probs.Count - 1;
    }
    
    public static double SampleDist(Random rnd, DataTypes.Distribution dist)
    {
        return dist switch
        {
            DataTypes.Distribution.Exponential => ExponentialDistribution.Sample(rnd, 1),
            DataTypes.Distribution.Normal => NormalDistribution.Sample(rnd, 1, 1),
            _ => 1
        };
    }
    
    public static BaseEventData GenerateCNEventData(Random rnd, Karyotype kar, CNEventP cnEventP)
    {
        List<(int id, long len)> seq = kar.ContigIds().Shuffle(rnd).Select(i => (i, kar.ContigLen(i))).ToList();

        switch (cnEventP.Type)
        {
            // Whole chromosome events
            case CNEventType.ChromDeletion:
            case CNEventType.ChromDuplication:
                return new ContigEventData(cnEventP, seq[0].id);
            
            case CNEventType.WholeGenomeDoubling:
                return new BaseEventData(cnEventP);
            
            // Tail events
            case CNEventType.TailDeletion:
            case CNEventType.BreakageFusionBridge:
                return new TailEventData(rnd, cnEventP, seq[0].id, seq[0].len);

            // Internal events
            case CNEventType.InternalDuplication:
            case CNEventType.InternalDeletion:
            case CNEventType.InternalInversion:
            case CNEventType.InvertedDuplication:
                return new InternalEventData(rnd, cnEventP, seq[0].id, seq[0].len);
            
            case CNEventType.Translocation:
                return new PairEventData(rnd, cnEventP, seq[0].id, seq[0].len, seq[1].id, seq[0].len);
            
            case CNEventType.Chromothripsis:
                return new ChromothripsisEventData(rnd, cnEventP, seq[0].id, seq[0].len);

            case CNEventType.Chromoplexy:
                return new ChromoplexyEventData(rnd, cnEventP, seq);

            case CNEventType.Pyrgo:
                return new PyrgoEventData(rnd, cnEventP, seq[0].id, seq[0].len);

            case CNEventType.Rigma:
                return new RigmaEventData(rnd, cnEventP, seq[0].id, seq[0].len);
            
            case CNEventType.TIChain:
            case CNEventType.TICycle:
            case CNEventType.TIBridge:
                return new TemplatedEventData(rnd, cnEventP, seq);

            default:
                throw new ArgumentOutOfRangeException(nameof(cnEventP.Type), cnEventP.Type, null);
        }
    }
}
