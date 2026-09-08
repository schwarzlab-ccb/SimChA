using SimChA.Data;
using SimChA.IO;

namespace SimChA.Computation;

public static class Fitness
{
    private const double EPSILON = 1e-8;

    private static double CalcTerm(double parameter, Func<double> termCalculation)
        => parameter > EPSILON ? termCalculation() * parameter : 0.0;

    public static double Calculate(Karyotype kar, RefGen refGen, FitParams fParams)
    {
        double stressTerm = CalcTerm(fParams.Stress, () => StressTerm(refGen.SexGenomeLen[(int) kar.Sex], kar.GenomeLen()));
        var geneData = refGen.SexGeneLists[(int)kar.Sex];
        double ogTerm = CalcTerm(fParams.TsgOg, () => TsgOgTerm(geneData[(int) GeneLT.OG], kar.GeneCounts[(int) GeneLT.OG]));
        double tsgTerm = CalcTerm(fParams.TsgOg, () => TsgOgTerm(geneData[(int) GeneLT.TSG], kar.GeneCounts[(int) GeneLT.TSG]));
        double essTerm = CalcTerm(fParams.Essentiality, () => EssTerm(geneData[(int) GeneLT.Ess], kar.GeneCounts[(int) GeneLT.Ess], fParams.HaploExponent));
        return 1 + stressTerm + ogTerm - tsgTerm + essTerm;
    }
    
    public static void LogCNs(IEnumerable<(Gene, int)> geneCNs)
    {
        Console.WriteLine("CNs:");
        foreach ((var gene, int cn) in geneCNs)
        {
            Console.WriteLine($"\tCN: {cn}; {gene}");
        }
    }

    public static double StressTerm(long refBaseCount, long baseCount)
        => Math.Min(0, 1 - baseCount / (double)refBaseCount);
    
    public static double TsgOgTerm(Gene[] genes, int[] geneCNs) 
        => genes.Select(gene=> Math.Log(1 + geneCNs[gene.GeneId])  / Math.Log(3) * gene.Score).Sum();

    public static double Zygosity(Gene[] genes, int[] geneCNs, int count)
        => genes.Sum(gene => geneCNs[gene.GeneId] == count ? 1 : 0);
    
    /// <summary>
    /// Copy number a gene's dosage is measured against. Fixed at the diploid 2 rather than read off
    /// the karyotype's current ploidy, which is the semantics of the CN==0 step this replaces: a
    /// WGD sample still has to lose all four copies before it is fully penalised. Making the
    /// reference ploidy-relative is a separate modelling decision, not part of this one.
    /// </summary>
    private const double DIPLOID = 2.0;

    /// <summary>
    /// The fraction of a gene's dosage that is missing, shaped by <paramref name="haploExponent"/>.
    /// CN >= 2 gives 0, CN == 1 gives 0.5^k and CN == 0 gives 1, so k alone sets what losing one of
    /// two copies costs relative to losing both:
    ///   k &lt; 1  one copy costs more than half of both -- strong haploinsufficiency
    ///   k == 1  linear in copies lost, one copy costs exactly half
    ///   k == 2  one copy costs a quarter (the default)
    ///   k -> inf  one copy costs nothing, recovering the CN==0-only step this replaces
    /// The last limit is what makes the change testable without touching code: a large exponent
    /// reproduces the previous model to within floating point.
    /// </summary>
    public static double DosageLoss(int geneCN, double haploExponent)
    {
        double lost = (DIPLOID - geneCN) / DIPLOID;
        return lost <= 0 ? 0.0 : Math.Pow(lost, haploExponent);
    }

    /// <summary>
    /// Essentiality penalty: negative, zero for a gene at its full diploid dosage.
    ///
    /// The step function this replaces (`min(CN-1, 0) * score`) was zero for every CN >= 1, so it
    /// only ever fired on outright homozygous deletion -- inactive in 89% of simulated samples,
    /// which is why Essentiality has been the least determined parameter of the fit in every run.
    /// </summary>
    public static double EssTerm(Gene[] genes, int[] geneCNs, double haploExponent)
        => -genes.Select(gene => DosageLoss(geneCNs[gene.GeneId], haploExponent) * gene.Score).Sum();

    /// <summary>
    /// Probability of accepting a proposed event, given the fitness it gains or loses and the
    /// acceptance offset delta. Glauber (Fermi) rather than the Metropolis `min(1, exp(dF - delta))`
    /// it replaces: the two agree wherever the proposal is clearly deleterious, but Metropolis
    /// clips to exactly 1 for every dF >= delta, and in that region the derivative with respect to
    /// delta -- and to every fitness parameter -- is exactly zero. At the September 2026 spice
    /// optimum (delta = 0.107) a third of accepted events sat in that clipped region, so a third of
    /// the run carried no information about the parameters being fitted. This form is strictly
    /// monotone in dF everywhere, so no event is ever selection-free.
    ///
    /// delta keeps its meaning as the fitness gain at which a proposal is accepted half the time,
    /// and the deleterious tail is unchanged: for dF - delta &lt;&lt; 0 both rules give exp(dF - delta).
    /// </summary>
    public static double AcceptProb(double deltaFitness, double acceptance)
        => 1.0 / (1.0 + Math.Exp(acceptance - deltaFitness));
}
