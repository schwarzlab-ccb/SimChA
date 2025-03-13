namespace SimChA.Data;

public class Gene(long start, long end, string chrom, string name, double deltaFitness, GeneLT listType)
    : GenRange(start, end, chrom)
{
    public string Name { get; } = name;
    public double DeltaFitness { get; } = deltaFitness;
    public GeneLT ListType { get; } = listType;
}