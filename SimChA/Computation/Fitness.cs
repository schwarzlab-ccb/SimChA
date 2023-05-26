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
               + FParams.TsgOg * (TsgOgTerm(ogCNs,karyotype.SexXX) - TsgOgTerm(tsgCNs,karyotype.SexXX)) 
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

    public static IEnumerable<(Gene, int)> CalcCNs(GeneListType geneListType, Karyotype karyotype)
    {
        var present = karyotype.GetPresentGenes(geneListType);
        var counts = present.GroupBy(g => g).ToDictionary(g =>g.Key, g => g.Count());
        var allSearched = GeneLists[geneListType].SelectMany(p => p.Value);
        var covered = allSearched.Where(g
            => !karyotype.IsMissing(g.Range) && (karyotype.SexXX || g.Range.ChrNo != ChrNo.chrY));
        return covered.Select(g => (g, counts.TryGetValue(g, out int count) ? count : 0));
    }
}
