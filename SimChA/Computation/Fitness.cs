// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Simulation;

namespace SimChA.Computation;

public static class Fitness
{
    public static double Calculate(
        Karyotype karyotype,
        Dictionary<ChrNo, List<Gene>> essentialGenes,
        Dictionary<ChrNo, List<Gene>> tsgOgGenes,
        SimParams simParams)
    {
        var tsgOgCNs = GeneCNs(karyotype.GetPresentGenes(tsgOgGenes));
        var essCNs = GeneCNs(karyotype.GetPresentGenes(essentialGenes));
        return 1 
               - simParams.StressFraction * StressTerm(karyotype.ContigCount) 
               + simParams.TsgOgFraction * TsgOgTerm(tsgOgCNs) 
               - simParams.EssentialFraction * EssTerm(essCNs);
    }

    // Represents the limitation of space in the nucleus - more chromosomes ==> more stress
    // TODO: This needs to be validated
    public static double StressTerm(int contigCount)
        => Math.Pow(Math.Max(0, contigCount - 46), 2);

    public static double EssTerm(IEnumerable<(Gene gene, int CN)> essCNs)
        => essCNs.Sum(g => Math.Max(1 - g.CN, 0) * g.gene.DeltaFitness);

    public static double TsgOgTerm(IEnumerable<(Gene gene, int CN)> tsgOgCNs)
        => tsgOgCNs.Sum(g => (g.CN - 2) * g.gene.DeltaFitness);

    private static IEnumerable<(Gene, int)> GeneCNs(IEnumerable<Gene> presentGenes)
        => presentGenes.GroupBy(g => g).Select(g => (g.Key, g.Count()));
}