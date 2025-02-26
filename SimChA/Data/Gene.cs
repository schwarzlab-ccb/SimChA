namespace SimChA.Data;

public class Gene : GenRange
{
    public string Name { get; }
    public double DeltaFitness { get; }
    
    public Gene(long start, long end, string chrom, string name, double deltaFitness) : base(start, end, chrom)
    {
        Name = name;
        DeltaFitness = deltaFitness;
    }
}