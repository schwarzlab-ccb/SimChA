// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.Computation;

public static class Fitness
{
    private static Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> GeneLists;
    private static FitnessParams FParams;

    public static double Calculate(Karyotype karyotype)
    {
        var tsgCNs = CalcCNs(GeneListType.TumorSuppressor, karyotype);
        var ogCNs = CalcCNs(GeneListType.Oncogene, karyotype);
        var essCNs = CalcCNs(GeneListType.Essentiality, karyotype);
        return 1 
               + FParams.Stress * StressTerm(karyotype.GenomeLen(), karyotype.SexXX) 
               + FParams.TsgOg * (TsgOgTerm(ogCNs) - TsgOgTerm(tsgCNs)) 
               + FParams.Essentiality * EssTerm(essCNs);
    }

    public static void SetStartingParams(Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> geneLists, FitnessParams fParams)
    {
        GeneLists = geneLists;
        FParams = fParams;
    }

    public static void LogCNs(IEnumerable<(Gene, int)> geneCNs)
    {
        Console.WriteLine("CNs:");
        foreach ((var gene, int cn) in geneCNs)
        {
            Console.WriteLine($"\tCN: {cn}; {gene}" );
        }
    }

    public static Dictionary<GeneListType, List<Gene>> GetGeneList(long start, long end, ChrNo chrNo)
    {
        if(GeneLists == null)
            return null;
        var geneList = new Dictionary<GeneListType, List<Gene>>();
        foreach(var gl in GeneLists.Keys)
        {
            geneList[gl] = GeneLists[gl][chrNo].FindAll(g => g.Range.IsInside(start, end, chrNo));
        }
        return geneList;
    }

    // Represents the limitation of space in the nucleus - more contigs ==> more stress
    // TODO: This needs to be validated
    public static double StressTerm(long baseCount, bool isFemale)
        => 1 - baseCount / (double) HGRef.GetGenomeLen(isFemale);

    public static double TsgOgTerm(IEnumerable<(Gene gene, int CN)> geneCNs)
        => geneCNs.Sum(g => (g.CN - 2) * g.gene.DeltaFitness);

    public static double EssTerm(IEnumerable<(Gene gene, int CN)> essCNs)
        => essCNs.Sum(g => Math.Min(g.CN - 1, 0) * g.gene.DeltaFitness);

    public static IEnumerable<(Gene, int)> CalcCNs(GeneListType geneListType, Karyotype karyotype)
    {
        var present = karyotype.GetPresentGenes(geneListType);
        var counts = present.GroupBy(g => g).ToDictionary(g =>g.Key, g => g.Count());
        var allSearched = GeneLists[geneListType].SelectMany(p => p.Value);
        var covered = allSearched.Where(g
            => !karyotype.IsMissing(g.Range) && (karyotype.SexXX || g.Range.ChrNo != ChrNo.chrY));
        return covered.Select(g => (g, counts.ContainsKey(g) ? counts[g] : 0));
    }
}