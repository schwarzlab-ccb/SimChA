using SimChA.Data;

namespace SimChA.Computation;

public static class CopyNumbers
{
    public static Dictionary<string, List<int>> GetJointSegmentation(List<string> chromNames, List<Sample> samples)
    {
        var breaks = samples.Select(s => s.Karyotype.CalcBreaks()).ToList();
        var segmentation = chromNames.ToDictionary(
            chrom => chrom, 
            chrom => breaks.SelectMany(br => br[chrom]).ToHashSet().OrderBy(val => val).ToList());
        return segmentation;
    }

    public static List<CopyNumber> CalcCNs(Karyotype karyotype, IDictionary<string, List<int>>? breaks = null) 
        => karyotype.CalcCNs(breaks ?? karyotype.CalcBreaks());

    public static double CalcPloidy(GenRef genRef, List<CopyNumber> copyNumbers, SexType sex)
        => 2 * copyNumbers.Select(c => c.Length * (c.CNH1 + c.CNH2)).Sum() 
           / (float) genRef.GetGenomeLen(sex);
}
