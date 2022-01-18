using SimChA.DataTypes;

namespace SimChA.Computation;

public static class SNPs
{
    public static List<SNP> CalculateSNPs(bool isFemale, int nrsnps, float hetrate = 1f)
    {
        List<SNP> snps = new();
        var referenceChromosomes = ReferenceGenome.GetChromosomes(isFemale);
        long totGenomeLength = ReferenceGenome.TotalLength(isFemale);
        Random random = new();

        foreach (var chrom in referenceChromosomes)
        {
            int curChromLength = ReferenceGenome.ChromosomeLengthMap[chrom];
            int nrChromSnps = (int) Math.Floor((double) nrsnps * curChromLength / totGenomeLength);
            for (int i = 0; i < nrChromSnps; i++)
            {
                snps.Add(new SNP(chrom, random.Next(1, curChromLength - 1), random.NextDouble() < hetrate));
            }
        }

        // Sort all SNPs
        snps = snps.OrderBy(snp => snp.AbsPos).ToList();
        return snps;
    }
}