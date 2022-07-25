using SimChA.DataTypes;

namespace SimChA.Computation;

public static class SNPBuilder
{
    public static List<SNP> CreateSNPs(Random random, bool isFemale, int nrsnps, float hetrate = 1f)
    {
        List<SNP> snps = new();
        var referenceChromosomes = ReferenceGenome.GetChromosomes(isFemale);
        long totGenomeLength = ReferenceGenome.TotalLength(isFemale);

        int snpId = 0;
        foreach (var chrom in referenceChromosomes)
        {
            int chromLength = ReferenceGenome.ChromosomeLengthMap[chrom];
            int nrChromSnps = (int) Math.Floor((double) nrsnps * chromLength / totGenomeLength);
            for (int i = 0; i < nrChromSnps; i++)
            {
                int pos = random.Next(1, chromLength - 1);
                bool isHet = random.NextDouble() < hetrate;
                snps.Add(new SNP(chrom, pos, isHet, snpId));
            }
            snpId += nrChromSnps;
        }

        snps = snps.OrderBy(snp => snp.AbsPos).ToList();
        return snps;
    }
}