using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.Computation;

public static class CopyNumbers
{
    public static List<CopyNumber> CalcCopyNumbers(Karyotype karyotype)
    {
        var allRegions = karyotype.GetAllRegions();
        var reference = ReferenceGenome.GetChromosomes(karyotype.IsFemale);
        var copyNumbers = CalcSegmentation(allRegions, reference);
        return copyNumbers;
        // TODO: Merge neighboring segments that have the same copy numbers
    }

    // Minimum consistent segmentation
    private static List<CopyNumber> CalcSegmentation(IList<Region> allRegs, IEnumerable<ChromNum> refChroms)
    {
        var result = new List<CopyNumber>();
        foreach (var refChrom in refChroms)
        {
            var curRegions = allRegs.Where(region => region.ChromId.ChromNum == refChrom).ToList();
            var curRegionsWithReference = curRegions.Append(ReferenceGenome.GetRegion(refChrom)).ToList();

            var starts = curRegionsWithReference.Select(r => r.Start);
            var ends = curRegionsWithReference.Select(r => r.End);
            var segmentBoundaries = starts.Concat(ends).Distinct().ToList();
            segmentBoundaries.Sort();

            var chromId = new ChromID(refChrom, true);
            for (int i = 0; i < segmentBoundaries.Count - 1; i++)
            {
                var seg = new Region(segmentBoundaries[i], segmentBoundaries[i + 1], chromId);
                var cn = new CopyNumber
                {
                    Segment = seg,
                    CNH1 = curRegions.Count(r => r.ChromId.Parent && seg.IsInside(r)),
                    CNH2 = curRegions.Count(r => !r.ChromId.Parent && seg.IsInside(r))
                };
                result.Add(cn);
            }
        }

        return result;
    }

    public static float CalcPloidy(List<CopyNumber> copyNumbers, bool isFemale)
    {
        long totalLength = ReferenceGenome.TotalLength(isFemale);
        float ploidy = copyNumbers
            .Select(c => (float)(c.Segment.End - c.Segment.Start) * (c.CNH1 + c.CNH2) / totalLength)
            .Sum();

        return ploidy;
    }

    private static string FirstLine(bool withSample, bool isFirst)
        => isFirst ? (withSample ? "sample_name\t" : "") + "chrom\tstart\tend\tcn_a\tcn_b\n" : "";

    public static string ToTSV(List<CopyNumber> copyNumbers, bool isFirst)
        => FirstLine(false, isFirst) + string.Join("\n", copyNumbers);

    public static string ToTSV(List<CopyNumber> copyNumbers, string sampleId, bool isFirst)
        => FirstLine(true, isFirst) + string.Join("\n", copyNumbers.Select(cn => $"{sampleId}\t{cn}"));
}