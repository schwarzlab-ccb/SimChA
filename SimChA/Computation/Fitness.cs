// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.Computation;

public static class Fitness
{
    public static double Calculate(
        Karyotype karyotype,
        GenRef genRef,
        FitnessParams fParams)
    {
        var tsgCNs = CalcCNs(genRef.GeneLists[GeneListType.TumorSuppressor], karyotype);
        var ogCNs = CalcCNs(genRef.GeneLists[GeneListType.Oncogene], karyotype);
        var essCNs = CalcCNs(genRef.GeneLists[GeneListType.Essentiality], karyotype);

        double stressTerm = StressTerm(genRef.GetGenomeLen(karyotype.SexXX), karyotype.GenomeLen())*fParams.Stress;
        double ogTerm = TsgOgTerm(genRef, ogCNs, karyotype.SexXX);
        double tsgTerm = TsgOgTerm(genRef, tsgCNs, karyotype.SexXX);
        double tsgogTerm = (ogTerm - tsgTerm)*fParams.TsgOg;
        double essTerm = EssTerm(genRef, essCNs, karyotype.SexXX, fParams.Haploinsufficiency)*fParams.Essentiality;
        
        return 1 + (stressTerm + tsgogTerm + essTerm)*fParams.TotalStrength;
    }

    public static double CalculateFromComponents(
        double stressTerm,
        double tsgogTerm,
        double essTerm,
        FitnessParams fParams)
    {
        return 1 + (stressTerm + tsgogTerm + essTerm)*fParams.TotalStrength;
    }

    public static void LogCNs(IEnumerable<(Gene, int)> geneCNs)
    {
        Console.WriteLine("CNs:");
        foreach ((var gene, int cn) in geneCNs)
        {
            Console.WriteLine($"\tCN: {cn}; {gene}" );
        }
    }

    public static double Sigmoid(double x)
        => 1 / (1 + Math.Exp(-((x * 1.5 - 1) * 10)));

    public static double Exponential(double x)
        => Math.Pow(x, 9) * 5;

    public static double Linear(double x)
        => x;

    public static double StressTerm(long refBaseCount, long baseCount)
        => Math.Min(0, 1 - baseCount / (double) refBaseCount);

    private static double ExpectedCN(GenRef genRef, string chrNo, bool sexXX)
    {
        if (chrNo == genRef.YChrName)
        {
            return sexXX ? 0 : 1;
        }
        if (chrNo == genRef.XChrName)
        {
            return sexXX ? 2 : 1;
        }
        return 2;
    }

    public static double TsgOgTerm(GenRef genRef, IEnumerable<(Gene gene, int CN)> geneCNs, bool sexXX)
    {
        var genesList = genRef.IncludeSexChromosomes ? geneCNs : geneCNs.Where(g => g.gene.Range.ChrNo != genRef.XChrName && g.gene.Range.ChrNo != genRef.YChrName);
        return genesList.Sum(g => (g.CN - ExpectedCN(genRef, g.gene.Range.ChrNo, sexXX)) * Linear(g.gene.DeltaFitness));
    }

    public static double EssTerm(GenRef genRef, IEnumerable<(Gene gene, int CN)> essCNs, bool sexXX, bool haploinsufficiency = false)
    {
        var genesList = genRef.IncludeSexChromosomes ? essCNs : essCNs.Where(g => g.gene.Range.ChrNo != genRef.XChrName && g.gene.Range.ChrNo != genRef.YChrName);
        return haploinsufficiency
            ? genesList.Sum(g => !(sexXX && g.gene.Range.ChrNo == genRef.YChrName) ? Math.Min(g.CN - ExpectedCN(genRef, g.gene.Range.ChrNo, sexXX), 0) * g.gene.DeltaFitness : 0)
            : genesList.Sum(g => !(sexXX && g.gene.Range.ChrNo == genRef.YChrName) ? Math.Min(g.CN - 1, 0) * g.gene.DeltaFitness : 0);
    }


    public static IEnumerable<(Gene, int)> CalcCNs(Dictionary<string, List<Gene>> searched, Karyotype karyotype)
    {
        var present = karyotype.GetPresentGenes(searched);
        var counts = present.GroupBy(g => g).ToDictionary(g =>g.Key, g => g.Count());
        var allSearched = searched.SelectMany(p => p.Value);
        return allSearched.Select(g => (g, counts.TryGetValue(g, out int count) ? count : 0));
    }
}
