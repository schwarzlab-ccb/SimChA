// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Simulation;

namespace SimChA.Computation;

public static class Fitness
{
    public static float Calculate(Karyotype karyotype, Dictionary<ChrNo, List<Gene>> essentialGenes,
        Dictionary<ChrNo, List<Gene>> tsgOgGenes, SimParams simParams)
    {
        float stress = CalcStress(karyotype.ContigCount);
        float tsgOgFitness = CalcTsgOgFitness(karyotype, tsgOgGenes);
        float essentialityFitness = CalcEssentialityFitness(karyotype, essentialGenes);
        return 1 - simParams.StressFraction * stress
            + simParams.TsgOgFraction * tsgOgFitness
            - simParams.EssentialFraction * essentialityFitness;
    }

    private static float CalcEssentialityFitness(Karyotype karyotype, Dictionary<ChrNo, List<Gene>> essentialGenes)
    {
        float essentialityFitness = 0;
        Dictionary<Gene, int> essentialOccurences =FindGeneMultiplications(essentialGenes, 
            karyotype.GetPresentGenes(essentialGenes));
        foreach(var essentialGene in essentialOccurences)
        {
            essentialityFitness += Math.Max(1 - essentialGene.Value, 0) * essentialGene.Key.DeltaFitness;
        }
        return essentialityFitness;   
    }

    private static float CalcTsgOgFitness(Karyotype karyotype, Dictionary<ChrNo, List<Gene>> tsgOgGenes)
    {
        float tsgOgFitness = 0;
        Dictionary<Gene, int> tsgOgOccurences = FindGeneMultiplications(tsgOgGenes, 
            karyotype.GetPresentGenes(tsgOgGenes));
        foreach(var tsgOgGenePresent in tsgOgOccurences)
        {
            tsgOgFitness += tsgOgGenePresent.Key.DeltaFitness * (tsgOgGenePresent.Value-2);
        }
        return tsgOgFitness;
    }

    private static Dictionary<Gene, int> FindGeneMultiplications(Dictionary<ChrNo, List<Gene>> geneLists, 
    List<Gene> presentGenes)
    {
        Dictionary<Gene, int> geneListMultiplication = new Dictionary<Gene, int>();
        foreach(var geneList in geneLists)
        {
            foreach(var gene in geneList.Value)
            {
                geneListMultiplication.Add(gene, presentGenes.Where(x => x == gene).Count());
            }
        }
        return geneListMultiplication;
    }

    // Represents the limitation of space in the nucleus - more chromosomes ==> more stress
    // TODO: This needs to be validated
    private static float CalcStress(int chrCount)
        => (float)Math.Pow(Math.Max(0, chrCount - 46), 2);
}