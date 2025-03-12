namespace SimChA.Data;

public class Gene(long start, long end, string chrom, string name, double deltaFitness)
    : GenRange(start, end, chrom)
{
    public string Name { get; } = name;
    public double DeltaFitness { get; } = deltaFitness;
}