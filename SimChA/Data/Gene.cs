namespace SimChA.Data;

public class Gene(long start, long end, string chrom, int listIndex, double deltaFitness, GeneLT listType)
    : GenRange(start, end, chrom)
{
    public double DeltaFitness { get; } = deltaFitness;
    public GeneLT ListType { get; } = listType;
}