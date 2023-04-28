using Extreme.Statistics.Distributions;
using SimChA.DataTypes;
using SimChA.Misc;

namespace SimChA.Simulation;

public static class Sampling
{
    public static long GetNormSeg(Random rnd, long contigLen, double meanFrac) 
        => Math.Max(1, Math.Min((long) Math.Round(contigLen * NormalDistribution.Sample(rnd, meanFrac, meanFrac / 3)), contigLen - 2));
    
    public static long GetExpSeg(Random rnd, long contigLen, long meanLen) 
        => Math.Min((long) Math.Round(ExponentialDistribution.Sample(rnd, meanLen)), contigLen - 2);
    
    public static long GetExpSeg(Random rnd, long contigLen, double meanFrac) 
        => Math.Min((long) Math.Round(contigLen * ExponentialDistribution.Sample(rnd, meanFrac)), contigLen - 2);

    // Get two positions within the contig (boundaries are excluded)
    public static long GetInternalPos(Random rnd, long contigLen)
        => rnd.NextInt64(1, contigLen - 1);
    
    public static int GetSiteCount(Random rnd, double frac)
        => (int) Math.Max(1, Math.Round(NormalDistribution.Sample(rnd, 1 / frac)));

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
        => Enumerable.Range(0, shardCount)
            .Select(_ => Sampling.GetInternalPos(rnd, contigLen))
            .Distinct()
            .OrderBy(i => i)
            .ToList();
    
    // TODO: Unify with internal operations
    public static Region CreateRandomRegion(Random rnd, List<Region> regions)
    {
        var region = regions.OrderBy(x => rnd.Next()).First();
        long start = rnd.NextInt64(region.Start, region.End - 1);
        long stop = rnd.NextInt64(start + 1, region.End);
        var newRegion = new Region(start, stop, region.ChrID, rnd.CoinFlip());
        return newRegion;
    }
    
    public static List<double> CreateRandomMixture(Random rnd, double[] concentrations)
        => concentrations.Any() ? new DirichletDistribution(concentrations).Sample(rnd).ToList() : new List<double>();

    
    public static CNEventP PickRandomEventP(Random rnd, List<CNEventP> eventPs)
    {
        var probs = eventPs.Select(ev => ev.Prob).ToList();
        probs = probs.Select(p => p / probs.Sum()).ToList();
        int index = Sampling.PickRandomIndex(rnd, probs);
        return eventPs[index];
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
}