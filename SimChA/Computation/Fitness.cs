// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.Data;
using SimChA.IO;

namespace SimChA.Computation;

public static class Fitness
{
    private const double EPSILON = 1e-8;

    public static double Calculate(
        Karyotype kar,
        GenRef genRef,
        FitParams fParams)
    {
        bool normGenes = fParams.GeneNormalization;
        double stressTerm = fParams.Stress > EPSILON
            ? StressTerm(genRef.GetGenomeLen(kar.Sex), kar.GenomeLen())
            : 0.0;
        double ogTerm = fParams.TsgOg > EPSILON 
            ? TsgOgTerm(genRef, CalcCNs(genRef.GeneLists[GeneListType.Oncogene], kar), kar.Sex, normGenes)
            : 0.0;
        double tsgTerm = fParams.TsgOg > EPSILON
            ? TsgOgTerm(genRef, CalcCNs(genRef.GeneLists[GeneListType.TumorSuppressor], kar), kar.Sex, normGenes)
            : 0.0;
        double essTerm = fParams.Essentiality > EPSILON
            ? EssTerm(genRef, CalcCNs(genRef.GeneLists[GeneListType.Essentiality], kar), kar.Sex, normGenes) 
            : 0.0;
        return 1 + stressTerm * fParams.Stress + (ogTerm - tsgTerm) * fParams.TsgOg + essTerm * fParams.Essentiality;
    }

    public static double CalculateFromComponents(
        double stressTerm,
        double tsgogTerm,
        double essTerm,
        FitParams fParams) 
        => 1 + stressTerm * fParams.Stress + tsgogTerm * fParams.TsgOg + essTerm * fParams.Essentiality;

    public static void LogCNs(IEnumerable<(Gene, int)> geneCNs)
    {
        Console.WriteLine("CNs:");
        foreach ((var gene, int cn) in geneCNs)
        {
            Console.WriteLine($"\tCN: {cn}; {gene}");
        }
    }

    public static double Sigmoid(double x)
        => 1 / (1 + Math.Exp(-((x * 1.5 - 1) * 10)));

    public static double Exponential(double x)
        => Math.Pow(x, 9) * 5;

    public static double Linear(double x)
        => x;

    public static double Tanh(double x)
        => Math.Tanh(x);

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
        => x < 0
            ? -Math.Log(1.0 - x / 2.0)
            : Math.Log(1.0 + x);

    // 0/0 => 1, i.e. the genes not present in the given sex contributed their default score
    private static double GetExpRatio(GenRef genRef, SexType sex, Gene g, int cn)
        => cn + (2 - ExpectedCN(genRef, g.Range.Chrom, sex));
    
    public static double TsgOgTerm(GenRef genRef, IEnumerable<(Gene gene, int CN)> geneCNs, SexType sex, bool normalizeGenes = false)
    {
        Func<(Gene gene, int CN), double> calsGene = 
            pair => Math.Log(1 + GetExpRatio(genRef, sex, pair.gene, pair.CN)) * pair.gene.DeltaFitness;

        return normalizeGenes 
            ? geneCNs.Average(calsGene)
            : geneCNs.Sum(calsGene);
    }

    private static bool IsAutosome(GenRef genRef, Gene gene) 
        => gene.Range.Chrom != genRef.XChrName && gene.Range.Chrom != genRef.YChrName;
    
    // @CODY TODO: This only works for autosomes currently.
    public static double Zygosity(GenRef genRef, IEnumerable<(Gene gene, int CN)> geneCNs, int count, bool normalizeGenes = false) 
        => normalizeGenes ?
            geneCNs.Average(pair => IsAutosome(genRef, pair.gene) && pair.CN == count ? 1 : 0) :
            geneCNs.Count(pair => IsAutosome(genRef, pair.gene) && pair.CN == count);

    public static double EssTerm(GenRef genRef, IEnumerable<(Gene gene, int CN)> essCNs, SexType sex, bool normalizeGenes = false)
    {
        var genesList = sex switch
        {
            SexType.Female => essCNs.Where(g => g.gene.Range.Chrom != genRef.YChrName),
            SexType.Male => essCNs,
            _ => essCNs.Where(pair => IsAutosome(genRef, pair.gene))
        };

        Func<(Gene gene, int CN), double> calsGene =
            pair => Math.Min(pair.CN - 1, 0) * pair.gene.DeltaFitness;
        
        return normalizeGenes
            ? genesList.Average(calsGene)
            : genesList.Sum(calsGene);
    }

    public static IEnumerable<(Gene, int)> CalcCNs(Dictionary<string, List<Gene>> searched, Karyotype karyotype)
    {
        var present = karyotype.GetPresentGenes(searched);
        var counts = present.GroupBy(g => g).ToDictionary(g => g.Key, g => g.Count());
        var allSearched = searched.SelectMany(p => p.Value);
        return allSearched.Select(g => (g, counts.GetValueOrDefault(g, 0)));
    }
}