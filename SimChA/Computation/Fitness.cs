// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.Computation;

public static class Fitness
{
    public static double Calculate(
        Karyotype karyotype,
        Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> geneLists,
        FitnessParams fParams)
    {
        var tsgCNs = CalcCNs(geneLists[GeneListType.TumorSuppressor], karyotype);
        var ogCNs = CalcCNs(geneLists[GeneListType.Oncogene], karyotype);
        var essCNs = CalcCNs(geneLists[GeneListType.Essentiality], karyotype);
        return 1 
               + fParams.Stress * StressTerm(karyotype.GenomeLen(), karyotype.SexXX) 
               + fParams.TsgOg * (TsgOgTerm(ogCNs, karyotype.SexXX) - TsgOgTerm(tsgCNs, karyotype.SexXX)) 
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

    // Represents the limitation of space in the nucleus - more contigs ==> more stress
    // TODO: This needs to be validated
    public static double StressTerm(long baseCount, bool isFemale)
        => 1 - baseCount / (double) HGRef.GetGenomeLen(isFemale);

    private static double ExpectedCN(ChrNo chrNo, bool sexXX)
        => chrNo switch
        {
            ChrNo.chrY => sexXX ? 0 : 1,
            ChrNo.chrX => sexXX ? 2 : 1,
            _ => 2
        };
    
    public static double TsgOgTerm(IEnumerable<(Gene gene, int CN)> geneCNs, bool sexXX)
        => geneCNs.Sum(g => (g.CN - ExpectedCN(g.gene.Range.ChrNo, sexXX)) * g.gene.DeltaFitness);

    public static double EssTerm(IEnumerable<(Gene gene, int CN)> essCNs)
        => essCNs.Sum(g => Math.Min(g.CN - 1, 0) * g.gene.DeltaFitness);

    public static IEnumerable<(Gene, int)> CalcCNs(Dictionary<ChrNo, List<Gene>> searched, Karyotype karyotype)
    {
        var present = karyotype.GetPresentGenes(searched);
        var counts = present.GroupBy(g => g).ToDictionary(g =>g.Key, g => g.Count());
        var allSearched = searched.SelectMany(p => p.Value);
        var covered = allSearched.Where(g => karyotype.SexXX || g.Range.ChrNo != ChrNo.chrY);
        return covered.Select(g => (g, counts.TryGetValue(g, out int count) ? count : 0));
    }
}