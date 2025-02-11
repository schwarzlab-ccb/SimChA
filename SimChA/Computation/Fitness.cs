// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.Data;
using SimChA.IO;

namespace SimChA.Computation;

public static class Fitness
{
    private const double EPSILON = 1e-8;

    public static double Calculate(
        Karyotype karyotype,
        GenRef genRef,
        FitParams fParams)
    {
        var tsgCNs = fParams.TsgOg > EPSILON
            ? CalcCNs(genRef.GeneLists[GeneListType.TumorSuppressor], karyotype)
            : new List<(Gene, int)>();
        var ogCNs = fParams.TsgOg > EPSILON
            ? CalcCNs(genRef.GeneLists[GeneListType.Oncogene], karyotype)
            : new List<(Gene, int)>();
        var essCNs = fParams.Essentiality > EPSILON
            ? CalcCNs(genRef.GeneLists[GeneListType.Essentiality], karyotype)
            : new List<(Gene, int)>();

        double stressTerm = fParams.Stress > EPSILON
            ? StressTerm(genRef.GetGenomeLen(karyotype.Sex), karyotype.GenomeLen())
            : 0.0;
        double ogTerm = TsgOgTerm(genRef, ogCNs, karyotype.Sex, fParams.GeneNormalization);
        double tsgTerm = TsgOgTerm(genRef, tsgCNs, karyotype.Sex, fParams.GeneNormalization);
        double essTerm = EssTerm(genRef, essCNs, karyotype.Sex, fParams.GeneNormalization);
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
        if (chrNo == genRef.YChrName)
        {
            return sex switch
            {
                SexType.Male => 1,
                _ => 0
            };
        }

        if (chrNo == genRef.XChrName)
        {
            return sex switch
            {
                SexType.Female => 2,
                SexType.Male => 1,
                _ => 0
            };
        }

        return 2;
    }

    public static double CountFn(double x)
        => x < 0
            ? -Math.Log(1.0 - x / 2.0)
            : Math.Log(1.0 + x);

    // 0/0 => 1, i.e. the genes not present in the given sex contributed their default score
    private static double GetExpRatio(GenRef genRef, SexType sex, Gene g, int cn)
        => ExpectedCN(genRef, g.Range.ChrNo, sex) > 0 ? cn / ExpectedCN(genRef, g.Range.ChrNo, sex) : 1;
    
    public static double TsgOgTerm(GenRef genRef, List<(Gene gene, int CN)> geneCNs, SexType sex, bool normalizeGenes = false)
    {
        if (!geneCNs.Any())
        {
            return 0;
        }
        int norm = normalizeGenes ? geneCNs.Count : 1;
        double sum = geneCNs.Sum(pair 
            => Math.Log2(1+GetExpRatio(genRef, sex, pair.gene, pair.CN)) * pair.gene.DeltaFitness) / norm;
        return sum;
    }

    public static double Zygosity(GenRef genRef, List<(Gene gene, int CN)> geneCNs, int count)
    {
        // TODO: This only works for autosomes currently.
        if (!geneCNs.Any())
        {
            return 0;
        }

        var genesList =
            geneCNs.Where(g => g.gene.Range.ChrNo != genRef.XChrName && g.gene.Range.ChrNo != genRef.YChrName).ToList();
        int zygosityCount = genesList.Count(g => g.CN == count);
        return zygosityCount / (double)genesList.Count;
    }

    public static double EssTerm(GenRef genRef, List<(Gene gene, int CN)> essCNs, SexType sex, bool normalizeGenes = false)
    {
        if (!essCNs.Any())
        {
            return 0;
        }

        var genesList = (sex switch
        {
            SexType.Female => essCNs.Where(g => g.gene.Range.ChrNo != genRef.YChrName),
            SexType.Male => essCNs,
            _ => essCNs.Where(g => g.gene.Range.ChrNo != genRef.XChrName && g.gene.Range.ChrNo != genRef.YChrName)
        }).ToList();
        int norm = normalizeGenes ? genesList.Count : 1;
        return genesList.Sum(g => Math.Min(g.CN - 1, 0) * g.gene.DeltaFitness) / norm;
    }

    public static List<(Gene, int)> CalcCNs(Dictionary<string, List<Gene>> searched, Karyotype karyotype)
    {
        var present = karyotype.GetPresentGenes(searched);
        var counts = present.GroupBy(g => g).ToDictionary(g => g.Key, g => g.Count());
        var allSearched = searched.SelectMany(p => p.Value);
        return allSearched.Select(g => (g, counts.GetValueOrDefault(g, 0))).ToList();
    }
}