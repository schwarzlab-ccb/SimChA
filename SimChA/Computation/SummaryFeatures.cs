using SimChA.DataTypes;
using SimChA.Simulation;
namespace SimChA.Computation;

public static class SummaryFeatures
{
    // Get the copy-number profiles of the clones in the sample
    public static Dictionary<Karyotype, List<CopyNumber>> GetCopyNumberProfiles(GenRef genRef, IList<Karyotype> kars)
    {
        var cnProfiles = new Dictionary<Karyotype, List<CopyNumber>>();
        foreach (var kar in kars)
        {
            var cnList = CopyNumbers.CalcCopyNumbers(genRef, kar, kar.SexXX).ToList();
            cnProfiles.Add(kar, cnList);
        }
        return cnProfiles;
    }

    public static List<long> GetSegLengths(GenRef genRef, IList<Karyotype> kars, bool includeCNNormal = false, bool includeLOH = false, bool includeSexChromosomes = false)
    {
        var segLengths = new List<long>();
        foreach (var cnProfile in GetCopyNumberProfiles(genRef, kars))
        {
            var cnList = cnProfile.Value;
            // Drop the sex chromosomes
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
            segLengths.AddRange(cnList.Select(cn => cn.Segment.Length));
        }
        return segLengths;
    }

    public static List<int> GetChangepoints(GenRef genRef, IList<Karyotype> kars, bool includeCNNormal = false, bool includeLOH = false, bool includeSexChromosomes = false)
    {
        var changepointList = new List<int>();

        foreach (var cnProfile in GetCopyNumberProfiles(genRef, kars))
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

    public static List<long> GetBreakpointsPerChromosome(GenRef genRef, IList<Karyotype> kars, bool includeSexChromosomes = false)
    {
        var breakpoints = new List<long>();
        foreach (var cnProfile in GetCopyNumberProfiles(genRef, kars))
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