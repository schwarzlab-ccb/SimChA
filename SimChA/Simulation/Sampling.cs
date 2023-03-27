using Extreme.Statistics.Distributions;
using SimChA.DataTypes;
using SimChA.Misc;

namespace SimChA.Simulation;

public static class Sampling
{
    public static double GetFraction(Random rnd, double mean)
        => Math.Clamp(ExponentialDistribution.Sample(rnd, mean), 0, 1);
    
    // Segment is at most 2 bases shorter than contig
    public static long GetSegLength(Random rnd, long contigLen, double mean) 
        => Math.Min((int)Math.Round(GetFraction(rnd, mean) * contigLen), contigLen - 2);

    // Get two positions within the contig (boundaries are excluded)
    public static long GetInternalPos(Random rnd, long contigLen)
        => rnd.NextInt64(1, contigLen - 1);
    
    // TODO: Needs better estimation
    public static int GetChromothripsisSiteCount(Random rnd, long contigLen)
        => rnd.Next(1, (int)Math.Pow(contigLen, 1 / 3f));
    
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
    public static Region CreateRandomRegion(List<Region> regions, Random rnd)
    {
        var region = regions.OrderBy(x => rnd.Next()).First();
        long start = rnd.NextInt64(region.Start, region.End - 1);
        long stop = rnd.NextInt64(start + 1, region.End);
        var newRegion = new Region(start, stop, region.ChrID, rnd.CoinFlip());
        return newRegion;
    }
}