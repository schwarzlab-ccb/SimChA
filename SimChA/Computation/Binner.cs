// Created by Dr. Cody B Duncan, 2023, cody.duncan@uk-koeln.de

using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.Computation;

public class Binner
{
    protected readonly Dictionary<string, List<long>> ChromosomeBins;
    protected readonly GenRef GenRef;
    public Binner(GenRef genRef)
    {
        GenRef = genRef;
        ChromosomeBins = GetChromosomeBins();
    }
    public Dictionary<string, List<CopyNumber>> GetBinnedCNProfiles(List<Sample> samples)
    {
        var binnedCNPs = new Dictionary<string, List<CopyNumber>>();
        foreach (var sample in samples)
        {
            foreach (var clone in sample.Clones)
            {
                binnedCNPs[sample.SampleId] = CopyNumbers.CalcBinnedCopyNumbers(GenRef, sample.Kars[clone.CloneId], ChromosomeBins, true).ToList();
            }
        }
        return binnedCNPs;
    }

    private Dictionary<string, List<long>> GetChromosomeBins()
    {
        var binSize = 1_000_000;
        var chromBins = new Dictionary<string, List<long>>();
        foreach (var chrom in GenRef.AllChrs)
        {
            var nFullBins = GenRef.ChrLengths[chrom] / binSize;
            var remainder = GenRef.ChrLengths[chrom] % binSize;
            // Adjusting the first and last bins
            var endBinSize = (long)(0.5 + remainder / 2.0);
            var offset = remainder - 2*endBinSize;
            var binList = new List<long>{0};
            for (int i = 0; i < nFullBins; i++)
            {
                binList.Add(i*binSize+endBinSize);
            }
            binList.Add(binList.Last()+endBinSize+offset-1);
            chromBins[chrom] = binList;
        }
        return chromBins;
    }
}