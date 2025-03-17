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
        double stressTerm = CalcTerm(fParams.Stress, () => StressTerm(genRef.GetGenomeLen(kar.Sex), kar.GenomeLen()));
        var geneData = genRef.GeneData[(int)kar.Sex];
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
    
    public static double TsgOgTerm(List<Gene> genes, List<int> geneCNs, bool normalizeGenes = false)
    {
        double sum = genes.Select((t, i) => Math.Log(1 + geneCNs[i]) * t.Score).Sum();
        return normalizeGenes ? sum / genes.Count : sum;
    }

    public static double Zygosity(List<int> geneCNs, int count, bool normalizeGenes = false)
    {
        double sum = geneCNs.Sum(i => i == count ? 1 : 0);
        return normalizeGenes ? sum / geneCNs.Count : sum;
    }
    
    public static double EssTerm(List<Gene> genes, List<int> geneCNs, bool normalizeGenes = false)
    {
        double sum = genes.Select((t, i) => Math.Min(geneCNs[i] - 1, 0) * t.Score).Sum();
        return normalizeGenes ? sum / genes.Count : sum;
    }
}
