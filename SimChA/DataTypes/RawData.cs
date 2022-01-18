using SimChA.Computation;

namespace SimChA.DataTypes;
using MathNet.Numerics.Distributions;

public class RawDataSingleSubclone
{
    private List<CopyNumber> _copynumbers;
    private float[] _baf;
    private float[] _logr;
    private float _purity;
    private float _ploidy;
    private float _baferror = 0.1F;
    private float _readstd = 2f;
    private int _readdepth = 10;
    private List<SNP> _snps;
    private int nrSNPs;

    public RawDataSingleSubclone(List<CopyNumber> copyNumbers, bool isFemale, float purity, float ploidy, int nrSNPsInitial = 100)
    {
        _copynumbers = copyNumbers;
        _purity = purity;
        _ploidy = ploidy;

        // Create SNPs
        _snps = SNPs.CalculateSNPs(isFemale, nrSNPsInitial);
        nrSNPs = _snps.Count;

        _baf = new float[nrSNPs];
        _logr = new float[nrSNPs];
    }

    public void CalcRawData()
    {
        int curSegmentId = 0;
        for (int i = 0; i < _snps.Count; i++)
        {
            curSegmentId = _copynumbers.FindIndex(curSegmentId, cn => 
                ReferenceGenome.ChromosomeAbsoluteStart(cn.Segment.ChromId.ChromNum + cn.Segment.End) < _snps[i].AbsPos);

            if (_copynumbers[curSegmentId].CNH1 + _copynumbers[curSegmentId].CNH2 <= 0) 
                continue;
            
            int readsH1 = (int)Math.Round(Math.Max(0, Normal.Sample(_copynumbers[curSegmentId].CNH1 * _readdepth, _readstd)));
            int readsH2 = (int)Math.Round(Math.Max(0, Normal.Sample(_copynumbers[curSegmentId].CNH2 * _readdepth, _readstd)));

            _logr[i] = (float)Math.Log2((readsH1 + readsH2) / (2f * _readdepth));
            _baf[i] = _snps[i].Heterozygous && readsH1+readsH2 > 0 ? (float) readsH2 / (readsH1 + readsH2) : -1;
            
            // TODO: Extra noise for logr and baf
        }
    }

    private string GetFirstLine(bool printHeader, string type)
        => printHeader ? $"\tchrom\tpos\t{type}\n" : "";
    
    public string BAFToTSV(bool printHeader = false)
    {
        var outputString = GetFirstLine(printHeader, "BAF");
        for (int i = 0; i < nrSNPs; i++) {
            outputString += $"snpID\t{_snps[i].Chrom}\t{_snps[i].Pos}\t{_baf[i]}\n";
        }
        return outputString;
    }
    public string logRToTSV(bool printHeader = false)
    {
        var outputString = GetFirstLine(printHeader, "logR");
        for (int i = 0; i < nrSNPs; i++) {
            outputString += $"snpID\t{_snps[i].Chrom}\t{_snps[i].Pos}\t{_logr[i]}\n";
        }
        return outputString;
    }
}