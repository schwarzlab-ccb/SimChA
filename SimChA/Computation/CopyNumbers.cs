using SimChA.DataTypes;

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
    private static List<CopyNumber> CalcSegmentation(IReadOnlyCollection<Region> allRegs, IEnumerable<ChromNum> refChroms)
    {
        var result = new List<CopyNumber>();
        foreach (var curRefChrom in refChroms)
        {
            var curRegions = allRegs.Where(region => region.ChromId.ChromNum == curRefChrom).ToList();
            curRegions.Add(ReferenceGenome.GetRegion(curRefChrom));

            var segmentBoundaries = curRegions
                .Select(r => r.Start)
                .Concat(curRegions.Select(r => r.End))
                .Distinct().ToList();
            segmentBoundaries.Sort();

            var chromId = new ChromID(curRefChrom, true);
            for (int i = 0; i < segmentBoundaries.Count - 1; i++)
            {
                var seg = new Region(segmentBoundaries[i], segmentBoundaries[i + 1], chromId);
                var cn = new CopyNumber
                {
                    Segment = seg,
                    // -1 because we add the reference region to curRegions above
                    CNH1 = curRegions.Count(r => r.ChromId.Parent && seg.IsInside(r)) - 1,
                    CNH2 = curRegions.Count(r => !r.ChromId.Parent && seg.IsInside(r))
                };
                result.Add(cn);
            }
        }
        return result;
    }

    private static string FirstLine(bool withSample, bool isFirst)
        => isFirst ? (withSample ? "sample_id\t" : "") + "chrom\tstart\tend\tcn_a\tcn_b\n" : "";

    public static string ToTSV(List<CopyNumber> copyNumbers, bool isFirst) 
        => FirstLine(false, isFirst) + string.Join("\n", copyNumbers);

    public static string ToTSV(List<CopyNumber> copyNumbers, string sampleId, bool isFirst) 
        => FirstLine(true, isFirst) + string.Join("\n", copyNumbers.Select(cn => $"{sampleId}\t{cn}"));
}