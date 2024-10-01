using SimChA.DataTypes;
using MathNet.Numerics.Statistics;
namespace SimChA.Computation;

public static class SummaryFeatures
{
    public static double GetSampleMajMinCN(List<CopyNumber> cnList, bool getMajor)
    {
        var weightedSum = cnList.Select(cn => cn.Segment.Length *
                                        (getMajor ? Math.Max(cn.CNH1, cn.CNH2)
                                                  : Math.Min(cn.CNH1, cn.CNH2)) ).ToList().Sum();
        var totalWeight = cnList.Select(cn => cn.Segment.Length).Sum();
        return weightedSum / (double)totalWeight;
    }
    public static (List<double> values, double max) GetMajMinCNs(Dictionary<string, List<CopyNumber>> cnProfiles, bool getMajor)
    {
        var cns = new List<double> ();
        var max = 0.0;
        foreach (var cnProfile in cnProfiles)
        {
            var cnList = cnProfile.Value;
            var mean = GetSampleMajMinCN(cnList, getMajor);
            if (mean > 10.0)
            {
                continue;
            }
            cns.Add(mean);
            if (mean > max) max = mean;
        }
        return (cns, max);
    }

    public static Dictionary<int, (double weight, List<double> segs)> GetStratifiedSegLengths(Dictionary<string, List<CopyNumber>> cnProfiles, 
        bool weightedByCount = false, bool includeCNNormal = true, bool includeLOH = true, bool includeSexChromosomes = false)
    {
        var cnLess  = new List<double>();
        var cnEqual = new List<double>();
        var cnMore  = new List<double>();
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
            cnLess.AddRange(cnList.Where(cn => cn.CNH1 + cn.CNH2 < 2).Select(cn => (double)cn.Segment.Length));
            cnEqual.AddRange(cnList.Where(cn => cn.CNH1 + cn.CNH2 == 2).Select(cn => (double)cn.Segment.Length));
            cnMore.AddRange(cnList.Where(cn => cn.CNH1 + cn.CNH2 > 2).Select(cn => (double)cn.Segment.Length));
        }
        var weights = new List<double>() {1.0/3, 1.0/3, 1.0/3};
        if (weightedByCount)
        {
            var nSegments = cnLess.Count + cnEqual.Count + cnMore.Count;
            weights = new List<double>() {cnLess.Count /(double) nSegments, cnEqual.Count /(double) nSegments, cnMore.Count /(double) nSegments};
        }
        var segLengths = new Dictionary<int, (double, List<double>)>
        {
            [0] = (weights[0], cnLess),
            [1] = (weights[1], cnEqual),
            [2] = (weights[2], cnMore)
        };
        return segLengths;
    }

    public static List<double> GetSegLengths(Dictionary<string, List<CopyNumber>> cnProfiles, 
        long cutoff = 20_000_000, bool includeCNNormal = false, bool includeLOH = false, bool includeSexChromosomes = false, bool weighted = false)
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
            if (weighted)
            {
                var lengths = cnList.Where(cn => cutoff <= 0 || cn.Segment.Length <= cutoff)
                                    .SelectMany(cn => Enumerable.Repeat((double)cn.Segment.Length, cn.CNH1 + cn.CNH2));
                segLengths.AddRange(lengths);
            }
            else
            {
                segLengths.AddRange(cnList.Where(cn => cutoff <= 0 || cn.Segment.Length <= cutoff)
                                          .Select(cn => (double) cn.Segment.Length));
            }
        }
        return segLengths;
    }

    public static List<double> GetMeanSegLength(Dictionary<string, List<CopyNumber>> cnps, bool includeCNNormal = false, bool includeLOH = false, bool includeSexChromosomes = false, bool weighted = false)
    {
        var meanSegLengths = new List<double>();
        foreach (var cnProfile in cnps)
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
            // Weighted means to account for copy-number of the segment
            if (weighted)
            {
                var lengths = cnList.SelectMany(cn => Enumerable.Repeat((double)cn.Segment.Length, cn.CNH1 + cn.CNH2));
                meanSegLengths.Add(lengths.DefaultIfEmpty(0).Average());
            }
            else
            {
                meanSegLengths.Add(cnList.Select(cn => (double) cn.Segment.Length).DefaultIfEmpty(0).Average());
            }
            
        }
        return meanSegLengths;
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

    public static Dictionary<string, List<int>> GetBreakpointsPerBin(GenRef genRef, Dictionary<string, List<CopyNumber>> cnProfiles, bool includeSexChromosomes, long SIZE)
    {
        var breakpoints = new Dictionary<string, List<int>>();
        var chrs = includeSexChromosomes ? genRef.AllChrs : genRef.ChrIDsForAutosomes();
        
        foreach (var cnProfile in cnProfiles)
        {
            var cnList = cnProfile.Value;
            var allBPs = new List<int>();
            foreach (var chr in chrs)
            {
                var chrLen = genRef.ChrLengths[chr];
                var chrSegs = cnList.Where(cn => cn.Segment.ChrNo == chr).ToList();
                var intervals = Enumerable.Range(0, (int)((chrLen + SIZE) / (double)SIZE) );
                var grouped = chrSegs.Where(cn => cn.Segment.End != chrLen) // Exclude endpoints of the chromosome
                                     .GroupBy(cn => (int)(cn.Segment.End / SIZE))
                                     .ToDictionary(g => g.Key, g => g.Count());

                var binned = Enumerable.Range(0, intervals.Count())
                                       .Select(i => grouped.ContainsKey(i) ? grouped[i] : 0)
                                       .ToList();
                allBPs.AddRange(binned);
            }
            breakpoints.Add(cnProfile.Key, allBPs);
        }
        return breakpoints;
    }

    public static List<double> GetBreakpointsDistribution(GenRef genRef, Dictionary<string, List<CopyNumber>> cnProfiles, bool includeSexChromosomes, bool perChrom = true, long SIZE = 10_000_000)
    {
        var bps = perChrom 
                    ? GetBreakpointsPerChromosome(genRef, cnProfiles, includeSexChromosomes)
                    : GetBreakpointsPerBin(genRef, cnProfiles, includeSexChromosomes, SIZE);
        // Assuming all vectors of the breakpoints for each sample have the same length (they should be)
        int listSize = bps.First().Value.Count;
        List<double> averages = Enumerable.Range(0, listSize)
            .Select(i => bps.Values.Select(list => list[i]).Average())
            .ToList();
        return averages;
    }

    public static double GetMeanBreakpoints(Dictionary<string, List<CopyNumber>> cnProfiles, Dictionary<string, int> eventCount, bool includeSexChromosomes)
        => GetBreakpoints(cnProfiles, eventCount, includeSexChromosomes).DefaultIfEmpty(0).Average();

    public static List<double> GetBreakpoints(Dictionary<string, List<CopyNumber>> cnProfiles, Dictionary<string, int> eventCount, bool includeSexChromosomes)
    {
        var breakpoints = new List<double>();
        foreach (var cnProfile in cnProfiles)
        {
            var cnList = cnProfile.Value;
            if (!includeSexChromosomes)
            {
                cnList = cnList.Where(cn => !(cn.Segment.ChrNo == "chrX" || cn.Segment.ChrNo == "chrY")).ToList();
            }
            double count = cnList.Count(cn => cn.Segment.End != cn.Segment.Length);
            if (eventCount[cnProfile.Key] > 0)
            {
                count /= eventCount[cnProfile.Key];
            }
            else
            {
                count = 0;
            }
            breakpoints.Add(count);
        }
        return breakpoints;
    }

    public static Dictionary<string, List<int>> GetBreakpointsPerChromosome(GenRef genRef, Dictionary<string, List<CopyNumber>> cnProfiles, bool includeSexChromosomes)
    {
        var breakpoints = new Dictionary<string, List<int>>();
        var chrs = includeSexChromosomes ? genRef.AllChrs : genRef.ChrIDsForAutosomes();
        foreach (var cnProfile in cnProfiles)
        {
            var cnList = cnProfile.Value;
            var allBPs = new List<int>();
            foreach (var chr in chrs)
            {
                var chrLen = genRef.ChrLengths[chr];
                var chrSegs = cnList.Where(cn => cn.Segment.ChrNo == chr).ToList();
                // Exclude endpoint of the chromosome
                var bps = chrSegs.Count(cn => cn.Segment.End != chrLen);
                allBPs.Add(bps);
            }
            breakpoints.Add(cnProfile.Key, allBPs);
        }
        return breakpoints;
    }

    public static List<double> GetHomozygousDeletionFraction(GenRef genRef, Dictionary<string, List<CopyNumber>> cnProfiles)
    {
        var fraction = new List<double>();
        foreach (var cnProfile in cnProfiles)
        {
            var homozygDelList = cnProfile.Value.Where(cn => cn.CNH1 + cn.CNH2 == 0 
                                                            && cn.Segment.ChrNo != "chrX" 
                                                            && cn.Segment.ChrNo != "chrY").ToList();
            var lenLost = homozygDelList.Select(cn => cn.Segment.Length).Sum();
            fraction.Add(lenLost/(double)genRef.AutosomeLen);
        }
        return fraction;
    }

    public static List<double> GetMeanCNAlongGenome(Dictionary<string, List<CopyNumber>> cnProfiles)
    {
        var meanCN = new List<double>();
        // Number of segments - assumes we already have binned copy-number profiles
        var cnLength = cnProfiles.First().Value.Count;

        for (int i = 0; i < cnLength; i++)
        {
            // The filtering process is for non-imputed CN profiles, since there may be bins with NaN values (e.g. missing segments)
            var filtered = cnProfiles.Where(kvp => kvp.Value[i].CNH1 + kvp.Value[i].CNH2 >= 0);
            var mean = filtered.Any() ? filtered.Average(kvp => kvp.Value[i].CNH1 + kvp.Value[i].CNH2) : 0.0;
            meanCN.Add(mean);
        }
        return meanCN;
    }

    public static double GetMeanPloidy(GenRef genRef, Dictionary<string, List<CopyNumber>> cnProfiles, Dictionary<string, SexEnum> sexes)
        => GetPloidy(genRef, cnProfiles, sexes, -1.0).DefaultIfEmpty(0).Average();

    public static List<double> GetPloidy(
        GenRef genRef, 
        Dictionary<string, List<CopyNumber>> cnProfiles,
        Dictionary<string, SexEnum> sexes, 
        double cutoff = 8.0)
        => cnProfiles.Select(kvp => CopyNumbers
                    .CalcPloidy(genRef, kvp.Value, sexes[kvp.Key]))
                    .Where(ploidy => cutoff <= 0 || ploidy <= cutoff)
                    .ToList();
                        
    // Produces a matrix of copy-number values for each chromosome in each sample
    public static Dictionary<string, Dictionary<string, double>> GetChrCopyNumberMatrix(
        List<string> chrs, 
        Dictionary<string, List<CopyNumber>> cnProfiles)
    {
        var matrix = new Dictionary<string, Dictionary<string, double>>();
        foreach (string chr in chrs)
        {
            var chrSpecificCN = new Dictionary<string, double>();
            foreach (var cnProfile in cnProfiles)
            {
                var chrCNPs = cnProfile.Value.Where(cn => cn.Segment.ChrNo == chr && cn.CNH1 + cn.CNH2 > 0).ToList();
                double val = 0.0;
                if (chrCNPs.Any())
                {
                    double weightedLen = chrCNPs.Select(cn => (cn.CNH1 + cn.CNH2)*cn.Segment.Length).Sum();
                    val = weightedLen / chrCNPs.Select(cn => cn.Segment.Length).Sum();
                }
                chrSpecificCN.Add(cnProfile.Key, val);
            }
            matrix.Add(chr, chrSpecificCN);
        }
        return matrix;
    }

    public static double GetMeanPloidy(Dictionary<string, Dictionary<string, double>> matrix)
        => matrix.SelectMany(row => row.Value)
                 .Select(col => col.Value)
                 .DefaultIfEmpty(0)
                 .Average();

    public static double GetMKV(Dictionary<string, Dictionary<string, double>> chrCNMatrix)
    {
        // Mean Karyotypic Variance is the variance of individual chromosomes across all samples
        // then averages across all chromosomes
        var chrVar = new Dictionary <string, double> ();
        foreach (string chr in chrCNMatrix.Keys)
        {
            var chrVals = chrCNMatrix[chr].Values;
            chrVar.Add(chr, chrVals.Variance());
        }
        return chrVar.Values.Average() / GetMeanPloidy(chrCNMatrix);
    }

    private static Dictionary<string, Dictionary<string, double>> TransposeChrCNMatrix( Dictionary<string, Dictionary<string, double>> chrCNMatrix)
    {
        var transposed = new Dictionary<string, Dictionary<string, double>>();
        foreach (var firstKey in chrCNMatrix.Keys)
        {
             foreach (var pair in chrCNMatrix[firstKey])
            {
                var secondKey = pair.Key;
                var value = pair.Value;
                if (!transposed.ContainsKey(secondKey))
                {
                    transposed[secondKey] = new Dictionary<string, double>();
                }
                transposed[secondKey][firstKey] = value;
            }
        }
        return transposed;
    }
    
    public static double GetAverageAneuploidy(Dictionary<string, Dictionary<string, double>> chrCNMatrix)
    {
        var matrix = TransposeChrCNMatrix(chrCNMatrix);
        var sampleVar = new Dictionary<string, double>();
        foreach (var sample in matrix.Keys)
        {
            var sampleVals = matrix[sample].Values;
            sampleVar.Add(sample, sampleVals.Variance());
        }
        return sampleVar.Values.Average();
    }
}