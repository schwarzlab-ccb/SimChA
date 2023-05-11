using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.Computation;

public static class CopyNumbers
{
    public static IEnumerable<CopyNumber> CalcCopyNumbers(Karyotype karyotype, bool isFemale)
    {
        var reference = HGRef.ChrIDsForSex(isFemale);
        var copyNumbers = reference.Select(c => 
            CalcChrCopyNumbers(karyotype.FindRegionsOfChr(c).ToList(), new ChrID(c, isFemale)));
        return copyNumbers.SelectMany(x => x).ToList();
        // TODO: Optimization: Merge neighboring segments that have the same copy numbers
    }

    // Calculate the segmentation of a chromosome based on the regions of the karyotype mapping to that chromosome
    private static IEnumerable<CopyNumber> CalcChrCopyNumbers(IReadOnlyCollection<Region> curRegs, ChrID id)
    {
        var result = new List<CopyNumber>();

        var curRegionsWithReference = curRegs.Append(HGRef.GetRegion(id.ChrNo)).ToList();

        var starts = curRegionsWithReference.Select(r => r.Start);
        var ends = curRegionsWithReference.Select(r => r.End);
        var segmentBoundaries = starts.Concat(ends).Distinct().ToList();
        segmentBoundaries.Sort();
        
        for (int i = 0; i < segmentBoundaries.Count - 1; i++)
        {
            var seg = new Region(segmentBoundaries[i], segmentBoundaries[i + 1], id);
            int cnh1 = curRegs.Count(r => r.ChrID.Parent && seg.IsInside(r));
            int cnh2 = curRegs.Count(r => !r.ChrID.Parent && seg.IsInside(r));
            var cn = new CopyNumber(seg, cnh1, cnh2);
            result.Add(cn);
        }

        return result;
    }

    public static double CalcPloidy(IEnumerable<CopyNumber> copyNumbers, bool isFemale)
    {
        long totalLength = HGRef.GetGenomeLen(isFemale) / 2;
        return copyNumbers.Select(c => (float)c.Segment.Length * (c.CNH1 + c.CNH2) / totalLength).Sum();
    }

    private static string Header(bool withSample, bool isFirst)
        => isFirst ? (withSample ? "sample_name\t" : "") + "chr\tstart\tend\tcn_a\tcn_b\n" : "";
    
    public static string ToTSV(IEnumerable<CopyNumber> copyNumbers, string sampleId, bool isFirst)
        => Header(true, isFirst) + string.Join("\n", copyNumbers.Select(cn => $"{sampleId}\t{cn.ToTSV()}"));
}