using SimChA.Data;

namespace SimChA.Computation;

public static class CopyNumbers
{
    public static IEnumerable<CopyNumber> CalcCNs(GenRef genRef, Karyotype karyotype) 
        => CalcCNs(genRef, karyotype, karyotype.CalcBreaks());

    public static IEnumerable<CopyNumber> CalcCNs(GenRef genRef, Karyotype karyotype, IDictionary<string, List<int>> breaks) 
        => genRef
            .ChrIDsForSex(karyotype.Sex)
            .SelectMany(c => karyotype.CalcChrCopyNumbers(breaks[c], c))
            .ToList();

    public static double CalcPloidy(GenRef genRef, IEnumerable<CopyNumber> copyNumbers, SexType sex)
    {
        long totalLength = genRef.GetGenomeLen(sex) / 2;
        return copyNumbers.Select(c => c.Length * (c.CNH1 + c.CNH2)).Sum() / (float) totalLength;
    }

    public static string Header(bool withSample)
        => (withSample ? "sample_name\t" : "") + "chr\tstart\tend\tcn_a\tcn_b\tn_snvs\n";
    
    public static string ToTSV(IEnumerable<CopyNumber> copyNumbers, string sampleId)
        => string.Join("\n", copyNumbers.Select(cn => $"{sampleId}\t{cn.ToTSV()}"));
}
