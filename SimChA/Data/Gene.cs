namespace SimChA.Data;

public class Gene(long start, long end, string chrom, string name, double deltaFitness, GeneListType listType)
    : GenRange(start, end, chrom)
{
    public string Name { get; } = name;
    public double DeltaFitness { get; } = deltaFitness;
    
    public GeneListType ListType { get; } = listType;
}