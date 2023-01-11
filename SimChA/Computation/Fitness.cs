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
        float stress = CalcStress(karyotype.ChrCount);
        Dictionary<Gene, int> tsgOgOccurences = FindeGeneMultiplications(tsgOgGenes, 
            karyotype.GetPresentGenes(tsgOgGenes));
        Dictionary<Gene, int> essentialOccurences = FindeGeneMultiplications(essentialGenes, 
            karyotype.GetPresentGenes(essentialGenes));
        float tsgOgFitness = 0;
        float essentialityFitness = 0;
        foreach(var tsgOgGenePresent in tsgOgOccurences)
        {
            tsgOgFitness += tsgOgGenePresent.Key.DeltaFitness * (tsgOgGenePresent.Value-2);
        }
        foreach(var essentialGene in essentialOccurences)
        {
            essentialityFitness = Math.Max(1 - essentialGene.Value, 0) * essentialGene.Key.DeltaFitness;
        }
        return 1 - simParams.StressFraction * stress
            + simParams.TsgOgFraction * tsgOgFitness
            - simParams.EssentialFraction * essentialityFitness;
    }

    private static Dictionary<Gene, int> FindeGeneMultiplications(Dictionary<ChrNo, List<Gene>> geneLists, 
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