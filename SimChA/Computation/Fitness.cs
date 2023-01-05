// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Simulation;

namespace SimChA.Computation;

public static class Fitness
{
    public static float Calculate(Karyotype karyotype, Dictionary<ChromNum, List<Gene>> essentialGenes,
        Dictionary<ChromNum, List<Gene>> tsgOgGenes, SimParams simParams)
    {
        float stress = CalcStress(simParams.StressFraction, karyotype.ChromCount);
        var essentialFound = karyotype.GetPresentGenes(essentialGenes);
        var tsgOgFound = karyotype.GetPresentGenes(tsgOgGenes);
        var essentialsMissing = FindMissingGenes(essentialFound, essentialGenes);
        var tsgOgMissing = FindMissingGenes(tsgOgFound, tsgOgGenes);
        var tsgOgCounts = tsgOgFound.GroupBy(x => x).ToList();
        float essentialityFitness = essentialsMissing.Sum(g => g.DeltaFitness);
        // Twice the value for missing genes (ploidy 0), -1 multiplicative factor for each missing gene (ploidy 1),
        // n - 2 for each overrepresented gene (ploidy 2+)
        float tsgOgFitness = 2 * tsgOgMissing.Sum(g => g.DeltaFitness) -
                             tsgOgCounts.Sum(g => g.Key.DeltaFitness * (g.Count() - 2));
        // parametrized linear combination of factors
        return stress + simParams.TsgOgFraction * tsgOgFitness + simParams.EssentialFraction * essentialityFitness;
    }

    private static IEnumerable<Gene> FindMissingGenes(IReadOnlyCollection<Gene> presentGenes,
        Dictionary<ChromNum, List<Gene>> geneList)
    {
        var missingGenes = new List<Gene>();
        foreach (var (_, allGenes) in geneList)
        {
            missingGenes.AddRange(allGenes.Except(presentGenes));
        }
        return missingGenes;
    }


    // Represents the limitation of space in the nucleus - more chromosomes ==> more stress
    // TODO: This needs to be validated
    private static float CalcStress(float stressFactor, int chromCount)
        => stressFactor * (float)Math.Pow(Math.Max(0, chromCount - 46), 2);
}