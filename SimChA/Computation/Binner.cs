// Created by Dr. Cody B Duncan, 2023, cody.duncan@uk-koeln.de

using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.Computation;

public class Binner
{
    protected readonly Dictionary<string, List<long>> ChromosomeBins;
    protected readonly GenRef GenRef;
    protected long BinWidth { get; }
    protected bool IncludeSexChromosomes { get; }
    public Binner(GenRef genRef, long binWidth = 1000000, bool includeSexChromosomes = false)
    {
        GenRef = genRef;
        BinWidth = binWidth;
        IncludeSexChromosomes = includeSexChromosomes;
        ChromosomeBins = GetChromosomeBins();
    }
    public Dictionary<string, List<CopyNumber>> GetBinnedCNProfiles(List<Sample> samples)
    {
        var binnedCNPs = new Dictionary<string, List<CopyNumber>>();
        foreach (var sample in samples)
        {
            foreach (var clone in sample.Clones)
            {
                binnedCNPs[sample.SampleId] = CopyNumbers.CalcBinnedCopyNumbers(sample.Kars[clone.CloneId], ChromosomeBins, true).ToList();
            }
        }
        return binnedCNPs;
    }

    public Dictionary<string, List<CopyNumber>> GetBinnedCNProfiles(Dictionary<string, Karyotype> karDict)
    {
        var binnedCNPs = new Dictionary<string, List<CopyNumber>>();
        foreach (var sample in karDict.Keys)
        {
            binnedCNPs[sample] = CopyNumbers.CalcBinnedCopyNumbers(karDict[sample], ChromosomeBins, true).ToList();
        }
        return binnedCNPs;
    }

    private Dictionary<string, List<long>> GetChromosomeBins()
    {
        var chromBins = new Dictionary<string, List<long>>();
        var chrs = IncludeSexChromosomes ? GenRef.AllChrs : GenRef.ChrIDsForAutosomes();
        foreach (var chrom in chrs)
        {
            var nFullBins = GenRef.ChrLengths[chrom] / BinWidth;
            var remainder = GenRef.ChrLengths[chrom] % BinWidth;
            // Adjusting the first and last bins
            var endBinSize = (long)(0.5 + remainder / 2.0);
            var offset = remainder - 2*endBinSize;
            var binList = new List<long>{0};
            for (int i = 0; i < nFullBins; i++)
            {
                binList.Add(i*BinWidth+endBinSize);
            }
            binList.Add(binList.Last()+endBinSize+offset-1);
            chromBins[chrom] = binList;
        }
        return chromBins;
    }
}