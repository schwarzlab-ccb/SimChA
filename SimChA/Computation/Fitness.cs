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
        return 1 
               + fParams.Stress * StressTerm(genRef.GetGenomeLen(karyotype.SexXX), karyotype.GenomeLen()) 
               + fParams.TsgOg * (TsgOgTerm(genRef, ogCNs, karyotype.SexXX) - TsgOgTerm(genRef, tsgCNs, karyotype.SexXX)) 
               + fParams.Essentiality * EssTerm(essCNs);
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

    // Represents the limitation of space in the nucleus - more contigs ==> more stress
    // TODO: This needs to be validated
    public static double StressTerm(long refBaseCount, long baseCount)
        => 1 - baseCount / (double) refBaseCount;

    private static double ExpectedCN(GenRef genRef, string chrNo, bool sexXX)
    {
        if (chrNo == genRef.XYChrName)
        {
            return sexXX ? 0 : 1;
        }
        if (chrNo == genRef.XXChrName)
        {
            return sexXX ? 2 : 1;
        }
        return 2;
    }

    public static double TsgOgTerm(GenRef genRef, IEnumerable<(Gene gene, int CN)> geneCNs, bool sexXX)
        => geneCNs.Sum(g => (g.CN - ExpectedCN(genRef, g.gene.Range.ChrNo, sexXX)) * Linear(g.gene.DeltaFitness));

    public static double EssTerm(IEnumerable<(Gene gene, int CN)> essCNs)
        => essCNs.Sum(g => Math.Min(g.CN - 1, 0) * g.gene.DeltaFitness);

    public static IEnumerable<(Gene, int)> CalcCNs(Dictionary<string, List<Gene>> searched, Karyotype karyotype)
    {
        var present = karyotype.GetPresentGenes(searched);
        var counts = present.GroupBy(g => g).ToDictionary(g =>g.Key, g => g.Count());
        var allSearched = searched.SelectMany(p => p.Value);
        return allSearched.Select(g => (g, counts.TryGetValue(g, out int count) ? count : 0));
    }
}