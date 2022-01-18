using SimChA.DataTypes;

namespace SimChA.Computation;

public struct SNPs
{
    public IEnumerable<ChromNum> ReferenceChromosomes;
    private List<SNP> _snps = new();
    private float _hetrate = 1F;
    private readonly Random _random = new();
    
    public SNPs(bool isFemale, int nrsnps)
    {
        ReferenceChromosomes = ReferenceGenome.GetChromosomes(isFemale);

        long totGenomeLength = ReferenceGenome.TotalLength(isFemale);

        foreach (var chrom in ReferenceChromosomes)
        {
            int curChromLength = ReferenceGenome.ChromosomeLengthMap[chrom];
            int nrChromSnps = (int) Math.Floor((double) nrsnps * curChromLength / totGenomeLength);
            for (int i = 0; i < nrChromSnps; i++)
            {
                _snps.Add(new SNP(chrom, _random.Next(1, curChromLength - 1), _random.NextDouble() < _hetrate));
            }
        }

        // Sort all SNPs
        _snps = _snps.OrderBy(snp => snp.AbsPos).ToList();
    }

    public override string ToString()
        => string.Join("\n", _snps.Select(r => r.ToString()));

    public List<SNP> AllSNPs() => _snps.ToList();

    public int Length => _snps.Count();
}