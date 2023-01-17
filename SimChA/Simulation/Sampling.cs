using Extreme.Statistics.Distributions;

namespace SimChA.Simulation;

public static class Sampling
{
    public static double GetFraction(Random rnd, double mean)
        => Math.Clamp(ExponentialDistribution.Sample(rnd, mean), 0, 1);
    
    // Segment is at most 2 bases shorter than contig
    public static int GetSegLength(Random rnd, int contigLen, double mean) 
        => Math.Min((int)Math.Round(GetFraction(rnd, mean) * contigLen), contigLen - 2);

    // Get two positions within the contig (boundaries are excluded)
    public static int GetInternalPos(Random rnd, int contigLen)
        => rnd.Next(1, contigLen - 1);
    
    // Needs better estimation
    public static int GetChromothripsisSiteCount(Random rnd, int contigLen)
        => rnd.Next(1, (int)Math.Pow(contigLen, 1 / 3f));
    
    // https://ashpublications.org/blood/article/134/Supplement_1/3767/424006/Chromoplexy-and-Chromothripsis-Are-Important
    private static int GetChromoplexySiteCount(Random rnd)
        => rnd.NextSingle() switch
        {
            var n when n < .46 => 3,
            var n when n < .64 => 4,
            var n when n < .74 => 5,
            var n when n < .79 => 6,
            _ => 2
        };

    public static List<int> GetStopsForShards(Random rnd, int contigLen, int shardCount) 
        => Enumerable.Range(0, shardCount)
            .Select(_ => Sampling.GetInternalPos(rnd, contigLen))
            .Distinct()
            .OrderBy(i => i)
            .ToList();
}