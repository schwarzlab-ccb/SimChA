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
        var geneRefCNs = refGen.SexGeneRefCNs[(int) kar.Sex];
        double essTerm = CalcTerm(fParams.Essentiality, () => EssTerm(geneData[(int) GeneLT.Ess],
            kar.GeneCounts[(int) GeneLT.Ess], geneRefCNs[(int) GeneLT.Ess]));
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
    /// Whether a karyotype has lost an essential gene outright.
    ///
    /// EvoSimulator rejects such a proposal rather than pricing it (see
    /// <see cref="IO.FitParams.ProhibitEssentialLoss"/>): losing both copies of a gene the cell
    /// cannot do without is inviability, not unfitness, and the empirical cohort shows it happening
    /// at a rate no penalty was reproducing -- 95 of 5.07 million sample-gene pairs across the 5958
    /// spice samples, 0.0019%, against 2.7x that in the best simulated cohort. Prohibiting it makes
    /// Essentiality mean one thing (the cost of haploinsufficiency) instead of setting the balance
    /// between two states at once.
    ///
    /// Short-circuits on the first hit; the scan is over the ~900 essential genes, whose counts
    /// Karyotype maintains incrementally.
    /// </summary>
    public static bool AnyEssentialLost(Gene[] genes, int[] geneCNs)
        => genes.Any(gene => geneCNs[gene.GeneId] == 0);

    /// <summary>
    /// Essentiality penalty: negative, and the score of every essential gene below diploid dosage.
    ///
    /// With <see cref="AnyEssentialLost"/> enforced this is a haploinsufficiency term -- CN == 0
    /// cannot arise, so the penalty is the scored count of genes down to a single copy, and
    /// FitParams.Essentiality is the price of exactly that. It replaces a graded
    /// `DosageLoss(CN)^k * score`, whose exponent set the CN==1 cost *relative to* CN==0; with
    /// CN==0 unreachable that ratio has no referent and k was exactly degenerate with Essentiality,
    /// so only their product was determined. Before that it was `min(CN-1, 0) * score`, which fired
    /// only on homozygous deletion and was inactive in 89% of simulated samples.
    ///
    /// CN == 0 is priced the same as CN == 1 rather than skipped, so the penalty stays monotone in
    /// copies lost for the modes that do not enforce the prohibition (MonteCarlo, FitnessMatching).
    /// Skipping it would make losing both copies cheaper than losing one.
    ///
    /// `refCNs` is the reference copy number per gene, parallel to `genes` -- two on the autosomes,
    /// one for a male's X and Y (see RefGen.SexGeneRefCNs). Scoring against a flat 2 charged every
    /// male for his own baseline on X.
    /// </summary>
    public static double EssTerm(Gene[] genes, int[] geneCNs, int[] refCNs)
        => -genes.Where((gene, idx) => geneCNs[gene.GeneId] < refCNs[idx]).Sum(gene => gene.Score);

    /// <summary>
    /// Probability of accepting a proposed event, given the fitness it gains or loses. Glauber
    /// (Fermi) rather than the Metropolis `min(1, exp(dF))` it replaces: the two agree wherever the
    /// proposal is clearly deleterious, but Metropolis clips to exactly 1 for every dF >= 0, and in
    /// that region the derivative with respect to every fitness parameter is exactly zero. At the
    /// September 2026 spice optimum a third of accepted events sat in that clipped region, so a
    /// third of the run carried no information about the parameters being fitted. This form is
    /// strictly monotone in dF everywhere, so no event is ever selection-free -- measured over the
    /// 2026-09-10 spice cohort, 99.3% of accepted events fall in the graded band
    /// 0.02 &lt; p &lt; 0.98 and only 0.6% saturate.
    ///
    /// There is no acceptance offset. `EvoParams.Acceptance` (delta) used to shift the midpoint and
    /// was removed because nothing could identify it. Two reasons, either sufficient: Stress, TsgOg
    /// and Essentiality enter dF linearly (see <see cref="Calculate"/>) and so already set the
    /// *scale* of the logistic, which leaves delta only a shift against three free scales; and with
    /// MaxTries at 100 a rejected proposal simply retries, so the acceptance rate reaches no
    /// observable at all -- the event count comes from the Gamma RateMean draw, not from how many
    /// proposals were refused. Measured over that cohort: 1.12 retries per accepted event on
    /// average, 17 at most, no slot anywhere near exhausting its budget, and moving delta from 0 to
    /// its last fitted value of 0.107 shifted mean acceptance 0.541 -> 0.517, a 4.5% uniform change
    /// with no differential signature. It had fitted with a CV of 0.37 and no correlate: a flat
    /// direction. Restoring it needs an observable first -- lowering MaxTries so that a refused
    /// proposal becomes a skip makes the realized event count depend on the acceptance rate, which
    /// is measurable.
    ///
    /// dF == 0 is therefore accepted exactly half the time: plain Glauber, the symmetric point.
    /// </summary>
    public static double AcceptProb(double deltaFitness)
        => 1.0 / (1.0 + Math.Exp(-deltaFitness));
}
