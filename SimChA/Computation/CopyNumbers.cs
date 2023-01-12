using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.Computation;

public static class CopyNumbers
{
    public static IEnumerable<CopyNumber> CalcCopyNumbers(Karyotype karyotype, bool isFemale)
    {
        var reference = ReferenceGenome.GetChromosomes(isFemale);
        var copyNumbers = reference.SelectMany(c => CalcSegmentation(karyotype.FindRegionsOfChr(c), c));
        return copyNumbers.ToList();
        // TODO: Optimization: Merge neighboring segments that have the same copy numbers
    }

    // Minimum consistent segmentation
    private static IEnumerable<CopyNumber> CalcSegmentation(IEnumerable<Region> curRegs, ChrNo refChrom)
    {
        var result = new List<CopyNumber>();

        var curRegionsWithReference = curRegs.Append(ReferenceGenome.GetRegion(refChrom)).ToList();

        var starts = curRegionsWithReference.Select(r => r.Start);
        var ends = curRegionsWithReference.Select(r => r.End);
        var segmentBoundaries = starts.Concat(ends).Distinct().ToList();
        segmentBoundaries.Sort();

        var chrID = new ChrID(refChrom, true);
        for (int i = 0; i < segmentBoundaries.Count - 1; i++)
        {
            var seg = new Region(segmentBoundaries[i], segmentBoundaries[i + 1], chrID);
            int cnh1 = curRegs.Count(r => r.ChrID.Parent && seg.IsInside(r));
            int cnh2 = curRegs.Count(r => !r.ChrID.Parent && seg.IsInside(r));
            var cn = new CopyNumber(seg, cnh1, cnh2);
            result.Add(cn);
        }


        return result;
    }

    public static float CalcPloidy(List<CopyNumber> copyNumbers, bool isFemale)
    {
        long totalLength = ReferenceGenome.TotalLength(isFemale);
        float ploidy = copyNumbers
            .Select(c => (float)c.Segment.Length * (c.CNH1 + c.CNH2) / totalLength)
            .Sum();

        return ploidy;
    }

    private static string Header(bool withSample, bool isFirst)
        => isFirst ? (withSample ? "sample_name\t" : "") + "chr\tstart\tend\tcn_a\tcn_b\n" : "";

    public static string ToTSV(IEnumerable<CopyNumber> copyNumbers, bool isFirst)
        => Header(false, isFirst) + string.Join("\n", copyNumbers.Select(cn => cn.ToTSV()));

    public static string ToTSV(IEnumerable<CopyNumber> copyNumbers, string sampleId, bool isFirst)
        => Header(true, isFirst) + string.Join("\n", copyNumbers.Select(cn => $"{sampleId}\t{cn.ToTSV()}"));
}