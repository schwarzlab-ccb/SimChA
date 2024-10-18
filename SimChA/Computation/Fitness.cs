// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.Computation;

public static class Fitness
{
    public static double Calculate(
        Karyotype karyotype,
        GenRef genRef,
        FitnessParams fParams)
    {
        var tsgCNs = CalcCNs(genRef.GeneLists[GeneListType.TumorSuppressor], karyotype);
        var ogCNs = CalcCNs(genRef.GeneLists[GeneListType.Oncogene], karyotype);
        var essCNs = CalcCNs(genRef.GeneLists[GeneListType.Essentiality], karyotype);

        double stressTerm = StressTerm(genRef.GetGenomeLen(karyotype.Sex), karyotype.GenomeLen())*fParams.Stress;
        double ogTerm = TsgOgTerm(genRef, ogCNs, karyotype.Sex, fParams.NormalizeGenes);
        double tsgTerm = TsgOgTerm(genRef, tsgCNs, karyotype.Sex, fParams.NormalizeGenes);
        double essTerm = EssTerm(genRef, essCNs, karyotype.Sex, fParams.NormalizeGenes, fParams.Haploinsufficiency)*fParams.Essentiality;
        double tsgogTerm = (ogTerm - tsgTerm)*fParams.TsgOg;
        
        
        return 1 + (stressTerm + tsgogTerm + essTerm)*fParams.TotalStrength;
    }

    public static double CalculateFromComponents(
        double stressTerm,
        double tsgogTerm,
        double essTerm,
        FitnessParams fParams)
    {
        return 1 + (stressTerm*fParams.Stress + tsgogTerm*fParams.TsgOg + essTerm*fParams.Essentiality)*fParams.TotalStrength;
    }

    public static void LogCNs(IEnumerable<(Gene, int)> geneCNs)
    {
        Console.WriteLine("CNs:");
        foreach ((var gene, int cn) in geneCNs)
        {
            Console.WriteLine($"\tCN: {cn}; {gene}" );
        }
    }

    public static double Sigmoid(double x)
        => 1 / (1 + Math.Exp(-((x * 1.5 - 1) * 10)));

    public static double Exponential(double x)
        => Math.Pow(x, 9) * 5;

    public static double Linear(double x)
        => x;

    public static double Tanh(double x)
        => Math.Tanh(x);

    public static double StressTerm(long refBaseCount, long baseCount)
        => Math.Min(0, 1 - baseCount / (double) refBaseCount);

    private static double ExpectedCN(GenRef genRef, string chrNo, SexEnum sex)
    {
        if (chrNo == genRef.YChrName)
        {
            return sex switch
            {
                SexEnum.Male => 1,
                _ => 0
            };
        }
        if (chrNo == genRef.XChrName)
        {
            return sex switch
            {
                SexEnum.Female => 2,
                SexEnum.Male => 1,
                _ => 0
            };
        }
        return 2;
    }

    /*public static double TsgOgTerm(GenRef genRef, IEnumerable<(Gene gene, int CN)> geneCNs, SexEnum sex, bool normalizeGenes = false)
    {
        var genesList = sex switch
        {
            SexEnum.Female => geneCNs.Where(g => g.gene.Range.ChrNo != genRef.YChrName),
            SexEnum.Male => geneCNs,
            _ => geneCNs.Where(g => g.gene.Range.ChrNo != genRef.XChrName && g.gene.Range.ChrNo != genRef.YChrName)
        };
        int norm = normalizeGenes ? genesList.Count() : 1;
        return genesList.Sum(g => (g.CN - ExpectedCN(genRef, g.gene.Range.ChrNo, sex)) * Linear(g.gene.DeltaFitness))/norm;
    }*/

    public static double TsgOgTerm(GenRef genRef, IEnumerable<(Gene gene, int CN)> geneCNs, SexEnum sex, bool normalizeGenes = false)
    {
        var genesList = sex switch
        {
            SexEnum.Female => geneCNs.Where(g => g.gene.Range.ChrNo != genRef.YChrName),
            SexEnum.Male => geneCNs,
            _ => geneCNs.Where(g => g.gene.Range.ChrNo != genRef.XChrName && g.gene.Range.ChrNo != genRef.YChrName)
        };
        int norm = normalizeGenes ? genesList.Count() : 1;

        return genesList.Sum(g => Math.Tanh(g.CN - ExpectedCN(genRef, g.gene.Range.ChrNo, sex)) * Linear(g.gene.DeltaFitness))/norm;
    }

    public static double OgTerm(GenRef genRef, IEnumerable<(Gene gene, int CN)> geneCNs, SexEnum sex, bool normalizeGenes = false)
    {
        var genesList = sex switch
        {
            SexEnum.Female => geneCNs.Where(g => g.gene.Range.ChrNo != genRef.YChrName),
            SexEnum.Male => geneCNs,
            _ => geneCNs.Where(g => g.gene.Range.ChrNo != genRef.XChrName && g.gene.Range.ChrNo != genRef.YChrName)
        };
        int norm = normalizeGenes ? genesList.Count() : 1;

        return genesList.Sum(g => (g.CN - ExpectedCN(genRef, g.gene.Range.ChrNo, sex)) * Linear(g.gene.DeltaFitness))/norm;
    }


    // TODO: Verify the sex here
    public static double EssTerm(GenRef genRef, IEnumerable<(Gene gene, int CN)> essCNs, SexEnum sex, bool normalizeGenes = false, bool haploinsufficiency = false)
    {
        var genesList = sex switch
        {
            SexEnum.Female => essCNs.Where(g => g.gene.Range.ChrNo != genRef.YChrName),
            SexEnum.Male => essCNs,
            _ => essCNs.Where(g => g.gene.Range.ChrNo != genRef.XChrName && g.gene.Range.ChrNo != genRef.YChrName)
        };
        int norm = normalizeGenes ? genesList.Count() : 1;
        return haploinsufficiency
            ? genesList.Sum(g => Math.Min(g.CN - ExpectedCN(genRef, g.gene.Range.ChrNo, sex), 0) * g.gene.DeltaFitness) / norm
            : genesList.Sum(g => Math.Min(g.CN - 1, 0) * g.gene.DeltaFitness) / norm;
    }


    public static IEnumerable<(Gene, int)> CalcCNs(Dictionary<string, List<Gene>> searched, Karyotype karyotype)
    {
        var present = karyotype.GetPresentGenes(searched);
        var counts = present.GroupBy(g => g).ToDictionary(g =>g.Key, g => g.Count());
        var allSearched = searched.SelectMany(p => p.Value);
        return allSearched.Select(g => (g, counts.TryGetValue(g, out int count) ? count : 0));
    }
}
