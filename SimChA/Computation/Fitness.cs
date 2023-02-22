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
               + fParams.Stress * StressTerm(karyotype.CountBases(), karyotype.SexXX) 
               + fParams.TsgOg * (TsgOgTerm(ogCNs) - TsgOgTerm(tsgCNs)) 
               + fParams.Essentiality * EssTerm(essCNs);
    }

    // Represents the limitation of space in the nucleus - more contigs ==> more stress
    // TODO: This needs to be validated
    public static double StressTerm(long baseCount, bool isFemale)
        => 1 - baseCount / (double) HGRef.GetGenomeLen(isFemale);

    public static double TsgOgTerm(IEnumerable<(Gene gene, int CN)> geneCNs)
        => geneCNs.Sum(g => (g.CN - 2) * g.gene.DeltaFitness);

    public static double EssTerm(IEnumerable<(Gene gene, int CN)> essCNs)
        => essCNs.Sum(g => Math.Min(g.CN - 1, 0) * g.gene.DeltaFitness);

    public static IEnumerable<(Gene, int)> CalcCNs(Dictionary<ChrNo, List<Gene>> searched, Karyotype karyotype)
    {
        var present = karyotype.GetPresentGenes(searched);
        var counts = present.GroupBy(g => g).ToDictionary(g =>g.Key, g => g.Count());
        var allSearched = searched.SelectMany(p => p.Value);
        return allSearched.Select(g => (g, counts.ContainsKey(g) ? counts[g] : 0));
    }
}