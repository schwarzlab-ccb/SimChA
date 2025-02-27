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
        double ogTerm = CalcTerm(fParams.TsgOg, () => TsgOgTerm(genRef, CalcCNs(genRef.GeneLists[GeneListType.Oncogene], kar), kar.Sex, normGenes));
        double tsgTerm = CalcTerm(fParams.TsgOg, () => TsgOgTerm(genRef, CalcCNs(genRef.GeneLists[GeneListType.TumorSuppressor], kar), kar.Sex, normGenes));
        double essTerm = CalcTerm(fParams.Essentiality, () => EssTerm(genRef, CalcCNs(genRef.GeneLists[GeneListType.Essentiality], kar), kar.Sex, normGenes));
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
    public static double CountFn(double x)
        => x < 0 ? -Math.Log(1.0 - x / 2.0) : Math.Log(1.0 + x);

    // 0/0 => 1, i.e. the genes not present in the given sex contributed their default score
    private static double GetExpRatio(GenRef genRef, SexType sex, Gene g, int cn)
        => cn + (2 - ExpectedCN(genRef, g.Chrom, sex));
    
    public static double TsgOgTerm(GenRef genRef, IEnumerable<(Gene gene, int CN)> geneCNs, SexType sex, bool normalizeGenes = false)
    {
        Func<(Gene gene, int CN), double> calsGene = 
            pair => Math.Log(1 + GetExpRatio(genRef, sex, pair.gene, pair.CN)) * pair.gene.DeltaFitness;

        return normalizeGenes 
            ? geneCNs.Average(calsGene)
            : geneCNs.Sum(calsGene);
    }

    private static bool IsAutosome(GenRef genRef, Gene gene) 
        => gene.Chrom != genRef.XChrName && gene.Chrom != genRef.YChrName;
    
    // @CODY TODO: This only works for autosomes currently.
    public static double Zygosity(GenRef genRef, IEnumerable<(Gene gene, int CN)> geneCNs, int count, bool normalizeGenes = false) 
        => normalizeGenes ?
            geneCNs.Average(pair => IsAutosome(genRef, pair.gene) && pair.CN == count ? 1 : 0) :
            geneCNs.Count(pair => IsAutosome(genRef, pair.gene) && pair.CN == count);

    private static IEnumerable<(Gene gene, int CN)> GenesForSex(GenRef genRef, IEnumerable<(Gene gene, int CN)> essCNs, SexType sex) 
        => sex switch
        {
            SexType.Female => essCNs.Where(g => g.gene.Chrom != genRef.YChrName),
            SexType.Male => essCNs,
            _ => essCNs.Where(pair => IsAutosome(genRef, pair.gene))
        };
    
    public static double EssTerm(GenRef genRef, IEnumerable<(Gene gene, int CN)> essCNs, SexType sex, bool normalizeGenes = false)
    {
        var genesList = GenesForSex(genRef, essCNs, sex);
        Func<(Gene gene, int CN), double> calsGene = pair => Math.Min(pair.CN - 1, 0) * pair.gene.DeltaFitness;
        return normalizeGenes ? genesList.Average(calsGene) : genesList.Sum(calsGene);
    }

    // TODO: The genes should exist on the karyotype actually
    public static IEnumerable<(Gene, int)> CalcCNs(Dictionary<string, List<Gene>> searched, Karyotype karyotype)
    {
        List<IEnumerable<(Gene, int)>> res = new();
        foreach ((string chrom, var geneList) in searched)
        {
            var present = karyotype.GetPresentGenes(chrom, geneList);
            var counts = present.ToLookup(g => g).ToDictionary(g => g.Key, g => g.Count());
            var allCounts = searched[chrom].Select(g => (g, counts.GetValueOrDefault(g.Name, 0)));
            res.Add(allCounts);
        }
        return res.SelectMany(v => v);
    }
}