using SimChA.Data;

namespace SimChA.Computation;

public static class CopyNumbers
{
    public static List<CopyNumber> CalcCNs(GenRef genRef, Karyotype karyotype)
    {
        karyotype.MergeRegions();
        return CalcCNs(genRef, karyotype, GetSegPoints(genRef, karyotype));
    }

    public static List<CopyNumber> CalcCNs(GenRef genRef, Karyotype karyotype, IDictionary<string, List<long>> breaks)
    {
        return genRef.ChrIDsForSex(karyotype.Sex)
            .SelectMany(c => CalcChrCopyNumbers(karyotype.FindRegionsOfChr(c).ToList(), breaks[c], c)).ToList();
    }
    
    private static IEnumerable<CopyNumber> CalcChrCopyNumbers(IReadOnlyCollection<Region> curRegs, List<long> breaks, string chrom)
    {
        var result = new List<CopyNumber>();
        for (int i = 0; i < breaks.Count - 1; i++)
        {
            long start = breaks[i];
            long end = breaks[i + 1];
            var seg = new GenRange(start, end, chrom);
            int cnh1 = curRegs.Count(r => r.Hap1 && seg.IsInsideOf(r));
            int cnh2 = curRegs.Count(r => !r.Hap1 && seg.IsInsideOf(r));
		    int nSNVs = curRegs.Sum(r => r.NumSNVsBetween(seg.Start, seg.End));
            var cn = new CopyNumber(seg, cnh1, cnh2, nSNVs);
            result.Add(cn);
        }
        return result;
    }

    public static double CalcPloidy(GenRef genRef, IEnumerable<CopyNumber> copyNumbers, SexType sex)
    {
        long totalLength = genRef.GetGenomeLen(sex) / 2;
        return copyNumbers.Select(c => c.Length * (c.CNH1 + c.CNH2)).Sum() / (float) totalLength;
    }

    public static string Header(bool withSample)
        => (withSample ? "sample_name\t" : "") + "chr\tstart\tend\tcn_a\tcn_b\tn_snvs\n";
    
    public static string ToTSV(IEnumerable<CopyNumber> copyNumbers, string sampleId)
        => string.Join("\n", copyNumbers.Select(cn => $"{sampleId}\t{cn.ToTSV()}"));

    private static List<long> GetSegPoints(GenRef genRef, string chrNo, Karyotype kar)
    {
        var segList = new HashSet<long> {0, genRef.ChrLengths[chrNo]};
        var segPoints = kar.FindRegionsOfChr(chrNo).SelectMany(r => new[] {r.AbsStart, r.AbsEnd});
        segList.UnionWith(segPoints);
        return segList.OrderBy(val => val).ToList();
    }

    public static Dictionary<string, List<long>> GetSegPoints(GenRef genRef, Karyotype kar) 
        => genRef.AllChrs.ToDictionary(chr => chr, chr => GetSegPoints(genRef, chr, kar));
}
