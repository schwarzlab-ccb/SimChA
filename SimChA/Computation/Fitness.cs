using SimChA.Data;
using SimChA.IO;

namespace SimChA.Computation;

public static class Fitness
{
    private const double EPSILON = 1e-8;

    private static double CalcTerm(double parameter, Func<double> termCalculation)
        => parameter > EPSILON ? termCalculation() * parameter : 0.0;

    public static double Calculate(Karyotype kar, GenRef genRef, FitParams fParams)
    {
        bool normGenes = fParams.GeneNormalization;
        double stressTerm = CalcTerm(fParams.Stress, () => StressTerm(genRef.SexGenomeLen[(int) kar.Sex], kar.GenomeLen()));
        var geneData = genRef.SexGeneLists[(int)kar.Sex];
        double ogTerm = CalcTerm(fParams.TsgOg, () => TsgOgTerm(geneData[(int) GeneLT.OG], kar.GeneCounts[(int) GeneLT.OG], normGenes));
        double tsgTerm = CalcTerm(fParams.TsgOg, () => TsgOgTerm(geneData[(int) GeneLT.TSG], kar.GeneCounts[(int) GeneLT.TSG], normGenes));
        double essTerm = CalcTerm(fParams.Essentiality, () => EssTerm(geneData[(int) GeneLT.Ess], kar.GeneCounts[(int) GeneLT.Ess], normGenes));
        return 1 + stressTerm + ogTerm - tsgTerm + essTerm;
    }
    
    public static void LogCNs(IEnumerable<(Gene, int)> geneCNs)
    {
        Console.WriteLine("CNs:");
        foreach ((var gene, int cn) in geneCNs)
        {
            Console.WriteLine($"\tCN: {cn}; {gene}");
        }
    }

    public static double StressTerm(long refBaseCount, long baseCount)
        => Math.Min(0, 1 - baseCount / (double)refBaseCount);
    
    public static double TsgOgTerm(Gene[] genes, int[] geneCNs, bool normalizeGenes = false)
    {
        double sum = genes.Select(gene=> Math.Log(1 + geneCNs[gene.GeneId]) * gene.Score).Sum();
        return normalizeGenes ? sum / genes.Length : sum;
    }

    public static double Zygosity(Gene[] genes, int[] geneCNs, int count, bool normalizeGenes = false)
    {
        double sum = genes.Sum(gene => geneCNs[gene.GeneId] == count ? 1 : 0);
        return normalizeGenes ? sum / genes.Length : sum;
    }
    
    public static double EssTerm(Gene[] genes, int[] geneCNs, bool normalizeGenes = false)
    {
        double sum = genes.Select(gene => Math.Min(geneCNs[gene.GeneId] - 1, 0) * gene.Score).Sum();
        return normalizeGenes ? sum / genes.Length : sum;
    }
}
