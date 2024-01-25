using SimChA.DataTypes;
using SimChA.Simulation;
using MathNet.Numerics.Statistics;
using SimChA.Optimization;
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
            cns.Add(mean);
            if (mean > max) max = mean;
        }
        return (cns, max);
    }
    public static (List<double> segs, double max) GetSegLengths(Dictionary<string, List<CopyNumber>> cnProfiles, bool includeCNNormal = false, bool includeLOH = false, bool includeSexChromosomes = false)
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
        var max = (segLengths.Count > 0) ? segLengths.Max() : 0;
        return (segLengths, max);
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

    public static (Dictionary<string, List<int>> values, int max) GetBreakpointsPerBin(GenRef genRef, Dictionary<string, List<CopyNumber>> cnProfiles, int SIZE = 10_000_000, bool includeSexChromosomes = false)
    {
        var breakpoints = new Dictionary<string, List<int>>();
        var maxBP = 0;
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
        return (breakpoints, maxBP);
    }

    public static List<double> GetBreakpointsDistribution(GenRef genRef, Dictionary<string, List<CopyNumber>> cnProfiles, bool includeSexChromosomes, int SIZE = 10_000_000)
    {
        var bps = GetBreakpointsPerChromosome(genRef, cnProfiles, includeSexChromosomes);
        // Assuming all vectors of the breakpoints for each sample have the same length (they should be)
        int listSize = bps.First().Value.Count;
        List<double> averages = Enumerable.Range(0, listSize)
            .Select(i => bps.Values.Select(list => list[i]).Average())
            .ToList();
        return averages;
    }

    public static Dictionary<string, List<int>> GetBreakpointsPerChromosome(GenRef genRef, Dictionary<string, List<CopyNumber>> cnProfiles, bool includeSexChromosomes = false)
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

    public static List<double> GetMeanCopyNumberAlongGenome(Dictionary<string, List<CopyNumber>> cnProfiles)
    {
        var meanCN = new List<double>();
        // Number of segments - assumes we already have binned copy-number profiles
        var cnLength = cnProfiles.First().Value.Count;

        for (int i = 0; i < cnLength; i++)
        {
            // TODO: What do we do with bins that are just empty?
            var filtered = cnProfiles.Where(kvp => kvp.Value[i].CNH1 + kvp.Value[i].CNH2 >= 0);
            var mean = filtered.Any() ? filtered.Average(kvp => kvp.Value[i].CNH1 + kvp.Value[i].CNH2) : 0.0;
            meanCN.Add(mean);
        }
        return meanCN;
    }

    public static (List<double> values, double max) GetPloidy(GenRef genRef, Dictionary<string, List<CopyNumber>> cnProfiles, Dictionary<string, bool> isFemaleDict, bool includeSexChromosomes = false)
    {
        var ploidies = includeSexChromosomes
                        ? cnProfiles.Select(kvp => CopyNumbers.CalcPloidy(genRef, kvp.Value, isFemaleDict[kvp.Key])).ToList()
                        : cnProfiles.Select(kvp => CopyNumbers.CalcAutosomePloidy(genRef, kvp.Value)).ToList();
        return (ploidies, ploidies.Max());
    }

    // Produces a matrix of copy-number values for each chromosome in each sample
    public static Dictionary<string, Dictionary<string, double>> GetChrCopyNumberMatrix(List<string> chrs, Dictionary<string, List<CopyNumber>> cnProfiles)
    {
        var matrix = new Dictionary<string, Dictionary<string, double>>();
        foreach (var chr in chrs)
        {
            var chrSpecificCN = new Dictionary<string, double>();
            foreach (var cnProfile in cnProfiles)
            {
                var chrCNPs = cnProfile.Value.Where(cn => cn.Segment.ChrNo == chr && cn.CNH1 + cn.CNH2 > 0).ToList();
                var val = 0.0;
                if (chrCNPs.Any())
                {
                    var weightedLen = chrCNPs.Select(cn => (cn.CNH1 + cn.CNH2)*cn.Segment.Length).Sum();
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
        foreach (var chr in chrCNMatrix.Keys)
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