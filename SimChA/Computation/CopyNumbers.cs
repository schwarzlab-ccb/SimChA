using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.Computation;

public static class CopyNumbers
{
    public static IEnumerable<CopyNumber> CalcCopyNumbers(GenRef genRef, Karyotype karyotype, bool isFemale)
    {
        var chrIDs = genRef.IncludeSexChromosomes ? genRef.ChrIDsForSex(isFemale) : genRef.ChrIDsForAutosomes();
        karyotype.MergeRegions();
        return chrIDs.SelectMany(c => CalcChrCopyNumbers(genRef, karyotype.FindRegionsOfChr(c), karyotype.GetMissingOfChr(c),c));
    } 
    
    public static IEnumerable<CopyNumber> CalcCopyNumbers(GenRef genRef, Karyotype karyotype, IDictionary<string, List<long>> segs, bool isFemale, bool keepMissing = false) 
    {
        var chrIDs = genRef.IncludeSexChromosomes ? genRef.ChrIDsForSex(isFemale) : genRef.ChrIDsForAutosomes();
        karyotype.MergeRegions();
        return chrIDs.SelectMany(c => CalcChrCopyNumbers(karyotype.FindRegionsOfChr(c).ToList(), karyotype.GetMissingOfChr(c), segs[c], c, keepMissing));
    }

    public static IEnumerable<CopyNumber> CalcBinnedCopyNumbers(Karyotype karyotype, IDictionary<string, List<long>> bins, bool keepMissing = false)
    {
        var chrIDs = bins.Keys;
        karyotype.MergeRegions();
        return chrIDs.SelectMany(c => CalcChrCopyNumbers(karyotype.FindRegionsOfChr(c).ToList(), karyotype.GetMissingOfChr(c), bins[c], c, keepMissing, false));
    }
    public static IEnumerable<CopyNumber> CalcChrCopyNumbers(GenRef genRef, IEnumerable<Region> curRegs, IList<GenRange> missing, string chrNo, bool keepMissing = false)
    {
        var regionList = curRegs.ToList();
        var starts = regionList.Select(r => r.Start).Append(0).Concat(missing.Select(r => r.Start));
        var ends = regionList.Select(r => r.End).Append(genRef.ChrLengths[chrNo]).Concat(missing.Select(r => r.End));
        var segmentBoundaries = starts.Concat(ends).Distinct().OrderBy(val => val).ToList();
        return CalcChrCopyNumbers(regionList, missing, segmentBoundaries, chrNo, keepMissing);
    }
    
    public static IEnumerable<CopyNumber> CalcChrCopyNumbers(IReadOnlyCollection<Region> curRegs, IList<GenRange> missing, IList<long> segs, string chrNo, bool keepMissing, bool joinSegments = true)
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
		        int nSNVs = curRegs.Sum(r => r.NumSNVsBetween(seg.Start, seg.End));
            	var cn = new CopyNumber(seg, cnh1, cnh2, nSNVs);
                result.Add(cn);
            }
            else if (keepMissing)
            {
                result.Add(new CopyNumber(seg, -1, -1, -1));
            }
        }
        if (joinSegments)
        {
            // Merge the result list such that adjacent segments with the same copy number are combined
            result = result.Aggregate(new List<CopyNumber>(), (acc, cn) =>
            {
                if (acc.Count == 0)
                {
                    acc.Add(cn);
                }
                else
                {
                    var last = acc[^1];
                    if (last.CNH1 == cn.CNH1 && last.CNH2 == cn.CNH2)
                    {
                        acc[^1] = new CopyNumber(new GenRange(last.Segment.Start, cn.Segment.End, last.Segment.ChrNo), last.CNH1, last.CNH2, last.NSNVs + cn.NSNVs);
                    }
                    else
                    {
                        acc.Add(cn);
                    }
                }
                return acc;
            });
        }
        return result;
    }

    public static double CalcPloidy(GenRef genRef, IEnumerable<CopyNumber> copyNumbers, bool isFemale)
    {
        long totalLength = genRef.GetGenomeLen(isFemale) / 2;
        return copyNumbers.Select(c => (float) c.Segment.Length * (c.CNH1 + c.CNH2) / totalLength).Sum();
    }

    public static double CalcAutosomePloidy(GenRef genRef, IEnumerable<CopyNumber> copyNumbers)
    {
        long totalLength = genRef.AutosomeLen / 2;
        return copyNumbers.Where(cn => cn.Segment.ChrNo != genRef.XChrName && cn.Segment.ChrNo != genRef.YChrName)
                          .Select(c => (float) c.Segment.Length * (c.CNH1 + c.CNH2) / totalLength).Sum();
    }

    private static string Header(bool withSample, bool isFirst)
        => isFirst ? (withSample ? "sample_name\t" : "") + "chr\tstart\tend\tcn_a\tcn_b\tn_snvs\n" : "";
    
    public static string ToTSV(IEnumerable<CopyNumber> copyNumbers, string sampleId, bool isFirst)
        => Header(true, isFirst) + string.Join("\n", copyNumbers.Select(cn => $"{sampleId}\t{cn.ToTSV()}"));

    public static List<long> GetSegPoints(GenRef genRef, string chrNo, IList<Karyotype> kars)
    {
        var segList = new HashSet<long> {0, genRef.ChrLengths[chrNo]};
        foreach (var kar in kars)
        {
            var segPoints = kar.FindRegionsOfChr(chrNo).SelectMany(r => new[] {r.Start, r.End});
            segList.UnionWith(segPoints);
        }
        return segList.OrderBy(val => val).ToList();
    }

    public static Dictionary<string, List<long>> GetSegPoints(GenRef genRef, IList<Karyotype> kars) 
        => genRef.AllChrs.ToDictionary(chr => chr, chr => GetSegPoints(genRef, chr, kars));
}
