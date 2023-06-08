using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.Computation;

public static class CopyNumbers
{
    public static IEnumerable<CopyNumber> CalcCopyNumbers(Karyotype karyotype, bool isFemale) 
        => HGRef.ChrIDsForSex(isFemale)
            .SelectMany(c => CalcChrCopyNumbers(karyotype.FindRegionsOfChr(c), karyotype.GetMissingOfChr(c),c));
    
    public static IEnumerable<CopyNumber> CalcCopyNumbers(Karyotype karyotype, IDictionary<ChrNo, List<long>> segs, bool isFemale) 
        => HGRef.ChrIDsForSex(isFemale)
            .SelectMany(c => CalcChrCopyNumbers(karyotype.FindRegionsOfChr(c).ToList(), karyotype.GetMissingOfChr(c), segs[c], c));

    public static IEnumerable<CopyNumber> CalcChrCopyNumbers(IEnumerable<Region> curRegs, IList<GenRange> missing, ChrNo chrNo)
    {
        var regionList = curRegs.ToList();
        var starts = regionList.Select(r => r.Start).Append(0).Concat(missing.Select(r => r.Start));
        var ends = regionList.Select(r => r.End).Append(HGRef.GetChromLen(chrNo)).Concat(missing.Select(r => r.End));
        var segmentBoundaries = starts.Concat(ends).Distinct().OrderBy(val => val).ToList();
        return CalcChrCopyNumbers(regionList, missing, segmentBoundaries, chrNo);;
    }
    
    public static IEnumerable<CopyNumber> CalcChrCopyNumbers(IReadOnlyCollection<Region> curRegs, IList<GenRange> missing, IList<long> segs, ChrNo chrNo)
    {
        var result = new List<CopyNumber>();
        for (int i = 0; i < segs.Count - 1; i++)
        {
            var seg = new GenRange(segs[i], segs[i + 1], chrNo);
            // Skip segments that are completely missing
            if (missing.All(m => !seg.IsInside(m)))
            {
                int cnh1 = curRegs.Count(r => r.Hap1 && seg.IsInside(r));
                int cnh2 = curRegs.Count(r => !r.Hap1 && seg.IsInside(r));
                var cn = new CopyNumber(seg, cnh1, cnh2);
                result.Add(cn);
            }
        }
        return result;
    }

    public static double CalcPloidy(IEnumerable<CopyNumber> copyNumbers, bool isFemale)
    {
        long totalLength = HGRef.GetGenomeLen(isFemale) / 2;
        return copyNumbers.Select(c => (float) c.Segment.Length * (c.CNH1 + c.CNH2) / totalLength).Sum();
    }

    private static string Header(bool withSample, bool isFirst)
        => isFirst ? (withSample ? "sample_name\t" : "") + "chr\tstart\tend\tcn_a\tcn_b\n" : "";

    public static string ToTSV(IEnumerable<CopyNumber> copyNumbers, string sampleId, bool isFirst)
        => Header(true, isFirst) + string.Join("\n", copyNumbers.Select(cn => $"{sampleId}\t{cn.ToTSV()}"));

    public static List<long> GetSegPoints(ChrNo chrNo, IEnumerable<Karyotype> kars)
    {
        var segList = new HashSet<long> {0, HGRef.GetChromLen(chrNo)};
        foreach (var kar in kars)
        {
            var segPoints = kar.FindRegionsOfChr(chrNo).SelectMany(r => new[] {r.Start, r.End});
            segList.UnionWith(segPoints);
        }
        return segList.OrderBy(val => val).ToList();
    }

    public static Dictionary<ChrNo, List<long>> GetSegPoints(IEnumerable<Karyotype> kars) 
        => Enum.GetValues<ChrNo>().ToDictionary(chr => chr, chr => GetSegPoints(chr, kars));
}