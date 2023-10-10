using SimChA.DataTypes;
using SimChA.Simulation;
using MathNet.Numerics.Statistics;
using SimChA.Optimization;
namespace SimChA.Computation;

public static class SummaryFeatures
{
    public static List<int> GetMinMajCNs(List<CopyNumber> cnList, bool isMajor)
        => cnList.Select(cn => isMajor ? Math.Max(cn.CNH1,cn.CNH2) : Math.Min(cn.CNH1, cn.CNH2)).ToList();

    public static (List<double> segList, long min, long max) GetSegLengthInfo(Dictionary<string, List<CopyNumber>> cnProfiles, bool includeCNNormal = false, bool includeLOH = false, bool includeSexChromosomes = false)
    {
        var segLengths = new List<double>();
        long minSeg = 1_000_000_000;
        long maxSeg = 0;
        foreach (var cnProfile in cnProfiles)
        {
            var cnList = cnProfile.Value;
            if (!includeSexChromosomes)
            {
                cnList = cnList.Where(cn => !(cn.Segment.ChrNo == "chrX" || cn.Segment.ChrNo == "chrY")).ToList();
            }
            if (!includeCNNormal)
            {
                cnList = cnList.Where(cn => !(cn.CNH1 == 1 && cn.CNH2 == 1)).ToList();
            }
            if (!includeLOH)
            {
                cnList = cnList.Where(cn => !(cn.CNH1 + cn.CNH2 == 2 && (cn.CNH1 == 0 || cn.CNH2 == 0))).ToList();
            }
            var minInList = cnList.Min(cn => cn.Segment.Length);
            if ( minInList < minSeg)
            {
                minSeg = minInList;
            }
            var maxInList = cnList.Max(cn => cn.Segment.Length);
            if (maxInList > maxSeg)
            {
                maxSeg = maxInList;
            }
            segLengths.AddRange(cnList.Select(cn => (double) cn.Segment.Length));
        }
        return (segLengths, minSeg, maxSeg);
    }

    public static List<int> GetChangepoints(Dictionary<string, List<CopyNumber>> cnProfiles, bool includeCNNormal = false, bool includeLOH = false, bool includeSexChromosomes = false)
    {
        var changepointList = new List<int>();

        foreach (var cnProfile in cnProfiles)
        {
            var cnList = cnProfile.Value;
            // Changepoint counts the step up or down between adjacent segments.
            // Left-most segment of copy-number profile uses a dummy diploid segment as its benchmark
            var leftSegmentCN = 2;
            var lastChr = null as string;
            foreach (var cn in cnList)
            {
                var thisChr = cn.Segment.ChrNo;
                if (thisChr != lastChr)
                {
                    leftSegmentCN = 2;
                    lastChr = thisChr;
                }
                if (!includeSexChromosomes && (thisChr == "chrX" || thisChr == "chrY"))
                {
                    continue;
                }
                var thisSegmentCN = cn.CNH1 + cn.CNH2;
                // By default, diploid segments are not counted
                if (thisSegmentCN == 2)
                {
                    if (!includeCNNormal && cn.CNH1 == 1 && cn.CNH2 == 1)
                    {
                        continue;
                    }
                    if (!includeLOH && (cn.CNH1 == 0 || cn.CNH2 == 0))
                    {
                        continue;
                    }
                }
                changepointList.Add(Math.Abs(leftSegmentCN - thisSegmentCN));
                leftSegmentCN = thisSegmentCN;
            }
        }
        return changepointList;
    }

    public static List<long> GetBreakpointsPerChromosome(Dictionary<string, List<CopyNumber>> cnProfiles, bool includeSexChromosomes = false)
    {
        var breakpoints = new List<long>();
        foreach (var cnProfile in cnProfiles)
        {
            var cnList = cnProfile.Value;
            var lastChr = null as string;
            var breakpointCount = 0;
            foreach (var cn in cnList)
            {
                var thisChr = cn.Segment.ChrNo;
                if (thisChr != lastChr)
                {
                    lastChr = thisChr;
                    if (breakpointCount > 0)
                    {
                        breakpoints.Add(breakpointCount);
                    }
                    breakpointCount = 0;
                }
                else
                {
                    breakpointCount += 1;
                }
            }
        }
        return breakpoints;
    }
}