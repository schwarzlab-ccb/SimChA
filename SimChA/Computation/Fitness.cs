// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Simulation;

namespace SimChA.Computation;

public static class Fitness
{
    public static double Calculate(
        Karyotype karyotype,
        Dictionary<ChrNo, List<Gene>> tsgGenes,
        Dictionary<ChrNo, List<Gene>> ogGenes, 
        Dictionary<ChrNo, List<Gene>> essentialGenes,
        SimParams simParams)
    {
        var tsgCNs = CalcCNs(tsgGenes, karyotype);
        var ogCNs = CalcCNs(ogGenes, karyotype);
        var essCNs = CalcCNs(essentialGenes, karyotype);
        return 1 
               - simParams.StressFraction * StressTerm(karyotype.ContigCount) 
               + simParams.TsgOgFraction * TsgOgTerm(tsgCNs)
               + simParams.TsgOgFraction * TsgOgTerm(ogCNs) 
               - simParams.EssentialFraction * EssTerm(essCNs);
    }

    // Represents the limitation of space in the nucleus - more contigs ==> more stress
    // TODO: This needs to be validated
    public static double StressTerm(int contigCount)
        => Math.Pow(Math.Max(0, contigCount - 46), 2);

    public static double TsgOgTerm(IEnumerable<(Gene gene, int CN)> geneCNs)
        => geneCNs.Sum(g => (g.CN - 2) * g.gene.DeltaFitness);

    public static double EssTerm(IEnumerable<(Gene gene, int CN)> essCNs)
        => essCNs.Sum(g => Math.Max(1 - g.CN, 0) * g.gene.DeltaFitness);

    public static IEnumerable<(Gene, int)> CalcCNs(Dictionary<ChrNo, List<Gene>> searched, Karyotype karyotype)
    {
        var present = karyotype.GetPresentGenes(searched);
        var counts = present.GroupBy(g => g).ToDictionary(g =>g.Key, g => g.Count());
        var allSearched = searched.SelectMany(p => p.Value);
        return allSearched.Select(g => (g, counts.ContainsKey(g) ? counts[g] : 0));
    }
}