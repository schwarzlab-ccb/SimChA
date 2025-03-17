namespace SimChA.Data;

public class Gene(long start, long end, string chrom, GeneLT listType, int geneId, double score)
    : GenRange(start, end, chrom)
{
    public double Score { get; } = score;
    public int GeneId { get; } = geneId;
    public GeneLT ListType { get; } = listType;
}