using SimChA.DataTypes;

namespace SimChA.Computation;

public static class SNPs
{
    public static List<SNP> CreateSNPs(Random random, bool isFemale, int nrsnps, float hetrate = 1f)
    {
        List<SNP> snps = new();
        var referenceChromosomes = ReferenceGenome.GetChromosomes(isFemale);
        long totGenomeLength = ReferenceGenome.TotalLength(isFemale);

        int curid = 0;
        foreach (var chrom in referenceChromosomes)
        {
            int curChromLength = ReferenceGenome.ChromosomeLengthMap[chrom];
            int nrChromSnps = (int) Math.Floor((double) nrsnps * curChromLength / totGenomeLength);
            for (int i = 0; i < nrChromSnps; i++)
            {
                snps.Add(new SNP(chrom, random.Next(1, curChromLength - 1), random.NextDouble() < hetrate, curid));
                curid++;
            }
        }

        snps = snps.OrderBy(snp => snp.AbsPos).ToList();
        return snps;
    }
}