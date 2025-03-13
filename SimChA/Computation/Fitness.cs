using SimChA.Data;
using SimChA.IO;

namespace SimChA.Computation;

public static class Fitness
{
    private const double EPSILON = 1e-8;

    private static double CalcTerm(double parameter, Func<double> termCalculation)
        => parameter > EPSILON ? termCalculation() * parameter : 0.0;

    public static double Calculate(Karyotype kar, GenRef genRef, FitParams fParams)
    {
        bool normGenes = fParams.GeneNormalization;
        double stressTerm = CalcTerm(fParams.Stress, () => StressTerm(genRef.GetGenomeLen(kar.Sex), kar.GenomeLen()));
        double ogTerm = CalcTerm(fParams.TsgOg, () => TsgOgTerm(genRef, kar.GeneCounts[GeneLT.OG], kar.Sex, normGenes));
        double tsgTerm = CalcTerm(fParams.TsgOg, () => TsgOgTerm(genRef, kar.GeneCounts[GeneLT.TSG], kar.Sex, normGenes));
        double essTerm = CalcTerm(fParams.Essentiality, () => EssTerm(kar.GeneCounts[GeneLT.Ess], normGenes));
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

    private static double ExpectedCN(GenRef genRef, string chrNo, SexType sex)
    {
        return (chrNo, sex) switch
        {
            ({ } chr, SexType.Male) when chr == genRef.YChrName || chr == genRef.XChrName  => 1,
            ({ } chr, SexType.Female) when chr == genRef.YChrName => 0,
            ({ } chr, SexType.Any) when chr == genRef.YChrName || chr == genRef.XChrName=> 0,
            _ => 2
        };
    }

    // 0/0 => 1, i.e. the genes not present in the given sex contributed their default score
    private static double GetExpRatio(GenRef genRef, SexType sex, Gene g, int cn)
        => cn + (2 - ExpectedCN(genRef, g.Chrom, sex));
    
    
    public static double TsgOgTerm(GenRef genRef, Dictionary<Gene, int> geneCNs, SexType sex, bool normalizeGenes = false)
    {
        Func<KeyValuePair<Gene, int>, double> calsGene = 
            pair => Math.Log(1 + GetExpRatio(genRef, sex, pair.Key, pair.Value)) * pair.Key.DeltaFitness;

        return normalizeGenes 
            ? geneCNs.Average(calsGene)
            : geneCNs.Sum(calsGene);
    }

    private static bool IsAutosome(GenRef genRef, Gene gene) 
        => gene.Chrom != genRef.XChrName && gene.Chrom != genRef.YChrName;
    
    public static double Zygosity(Dictionary<Gene, int> geneCNs, int count, bool normalizeGenes = false) 
        => normalizeGenes ?
            geneCNs.Average(pair => pair.Value == count ? 1 : 0) : geneCNs.Count(pair => pair.Value == count);
    
    public static double EssTerm(Dictionary<Gene, int> geneCNs, bool normalizeGenes = false)
    {
        Func<KeyValuePair<Gene, int>, double> calsGene = 
            pair => Math.Min(pair.Value - 1, 0) * pair.Key.DeltaFitness;
        return normalizeGenes ? geneCNs.Average(calsGene) : geneCNs.Sum(calsGene);
    }
}