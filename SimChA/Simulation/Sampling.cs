using Extreme.Statistics.Distributions;

namespace SimChA.Simulation;

public static class Sampling
{
    public static double GetFraction(Random rnd, double mean)
        => Math.Clamp(ExponentialDistribution.Sample(rnd, mean), 0, 1);
    
    // Segment is at most 2 bases shorter than contig
    public static int GetSegLength(double fraction, Contig contig) 
        => Math.Min((int)Math.Round(fraction * contig.Length()), contig.Length() - 2);
    
    // Get two positions within the contig (boundaries are excluded)
    public static (int start, int end) GetInternalRange(Random rnd, Contig contig, int segLength)
    {
        int start = rnd.Next(contig.Length() - segLength + 1);
        int end = Math.Min(start + segLength + 1, contig.Length() - 1);
        return (start, end);
    }
    
    // Get two positions within the contig (boundaries are excluded)
    public static int GetUniformPos(Random rnd, Contig contig)
        => rnd.Next(1, contig.Length() - 1);

    private static (int, bool) GetTail(Random rnd, int segLength, Contig contig, bool fiveToThree)
    {
        int pos = fiveToThree ? segLength - 1 : contig.Length() - segLength - 1;
        return (pos, fiveToThree);
    }
    
    // Needs better estimation
    public static int GetChromothripsisSiteCount(Random rnd, Contig contig)
        => rnd.Next(1, (int)Math.Pow(contig.Length(), 1 / 3f));
}