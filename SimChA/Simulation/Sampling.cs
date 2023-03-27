using Extreme.Statistics.Distributions;

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
    
    // Needs better estimation
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
    

    // Only two to seven regions were added therefore max of 7 -> https://doi.org/10.1038/s41586-019-1913-9
    public static int BetaDistribution(Random rnd)
    {
        double u1 = rnd.NextDouble();
        double u2 = rnd.NextDouble();

        double x = Math.Pow(u1, 1.0/0.1);
        double y = Math.Pow(u2, 1.0/7.0);

        return (int)Math.Round((x/(x+y)*5)+2);
    }

    public static long LongRandom(long min, long max, Random rnd)
    {
        byte[] buf = new byte[8];
        rnd.NextBytes(buf);
        long longRand = BitConverter.ToInt64(buf, 0);
        return (Math.Abs(longRand % (max-min)) + min);
    }
}