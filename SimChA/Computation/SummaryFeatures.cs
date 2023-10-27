using SimChA.DataTypes;
using SimChA.Simulation;
using MathNet.Numerics.Statistics;
using SimChA.Optimization;
namespace SimChA.Computation;

public static class SummaryFeatures
{
    public static double GetSampleMajMinCN(List<CopyNumber> cnList, bool getMajor)
        => cnList.Select(cn => getMajor 
                        ? Math.Max(cn.CNH1,cn.CNH2)
                        : Math.Min(cn.CNH1, cn.CNH2) ).ToList().Average();

    public static (List<double> values, double max) GetMajMinCNs(Dictionary<string, List<CopyNumber>> cnProfiles, bool getMajor)
    {
        var cns = new List<double> ();
        var max = 0.0;
        foreach (var cnProfile in cnProfiles)
        {
            var cnList = cnProfile.Value;
            var mean = GetSampleMajMinCN(cnList, getMajor);
            cns.Add(mean);
            if (mean > max) max = mean;
        }
        return (cns, max);
    }
    public static (List<double>, double max) GetSegLengths(Dictionary<string, List<CopyNumber>> cnProfiles, bool includeCNNormal = false, bool includeLOH = false, bool includeSexChromosomes = false)
    {
        var segLengths = new List<double>();
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
            segLengths.AddRange(cnList.Select(cn => (double) cn.Segment.Length));
        }
        return (segLengths, segLengths.Max());
    }

    public static (List<double> values, int max) GetChangepointInfo(Dictionary<string, List<CopyNumber>> cnProfiles, bool includeCNNormal = false, bool includeLOH = false, bool includeSexChromosomes = false)
    {
        var changepointList = new List<double>();
        var maxChange = 0;
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
                var change = Math.Abs(leftSegmentCN - thisSegmentCN);
                if (change > maxChange)
                {
                    maxChange = change;
                }
                changepointList.Add(change);
                leftSegmentCN = thisSegmentCN;
            }
        }
        return (changepointList, maxChange);
    }

    public static (List<double> values, int max) GetBreakpointsPerChromosome(Dictionary<string, List<CopyNumber>> cnProfiles, bool includeSexChromosomes = false)
    {
        var breakpoints = new List<double>();
        var maxBP = 0;
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
                    if (lastChr != null)
                    {
                        if (breakpointCount > maxBP)
                        {
                            maxBP = breakpointCount;
                        }
                        breakpoints.Add(breakpointCount);
                    }
                    lastChr = thisChr;
                    breakpointCount = 0;
                }
                else
                {
                    breakpointCount += 1;
                }
            }
        }
        return (breakpoints, maxBP);
    }

    public static (List<double> values, double max) GetHomozygousDeletionFraction(Dictionary<string, List<CopyNumber>> cnProfiles, long autosomeLength)
    {
        var fraction = new List<double>();
        foreach (var cnProfile in cnProfiles)
        {
            var homozygDelList = cnProfile.Value.Where(cn => cn.CNH1 + cn.CNH2 == 0 
                                                            && cn.Segment.ChrNo != "chrX" 
                                                            && cn.Segment.ChrNo != "chrY").ToList();
            fraction.AddRange(homozygDelList.Select(cn => (double) cn.Segment.Length/autosomeLength));
        }

        return (fraction, fraction.Max());
    }

    public static (List<double> values, double max) GetMeanCopyNumberAlongGenome(Dictionary<string, List<CopyNumber>> cnProfile)
    {
        var meanCN = new List<double>();
        var max = 0;

        return (meanCN, max);
    }
}