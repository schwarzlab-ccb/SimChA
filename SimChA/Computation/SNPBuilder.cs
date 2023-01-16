using SimChA.DataTypes;

namespace SimChA.Computation;

public static class SNPBuilder
{
    public static List<SNP> CreateSNPs(Random random, bool isFemale, int nrsnps, float hetrate = 1f)
    {
        List<SNP> snps = new();
        var referenceChrs = ReferenceGenome.GetChromosomes(isFemale);
        long totGenomeLength = ReferenceGenome.TotalLength(isFemale) / 2;

        int snpId = 0;
        foreach (var chr in referenceChrs)
        {
            int chrLength = ReferenceGenome.ChromosomeLengthMap[chr];
            int nrChrSnps = (int)Math.Floor((double)nrsnps * chrLength / totGenomeLength);
            for (int i = 0; i < nrChrSnps; i++)
            {
                int pos = random.Next(1, chrLength - 1);
                bool isHet = random.NextDouble() < hetrate;
                snps.Add(new SNP(chr, pos, isHet, snpId));
            }

            snpId += nrChrSnps;
        }

        snps = snps.OrderBy(snp => snp.AbsPos).ToList();
        return snps;
    }
}