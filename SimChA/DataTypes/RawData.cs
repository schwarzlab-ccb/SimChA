namespace SimChA.DataTypes;
using MathNet.Numerics.Distributions;

public class RawDataSingleSubclone
{
    private CopyNumbers _copynumbers;
    private float[] _baf;
    private float[] _logr;
    private float _purity;
    private float _ploidy;
    private float _baferror = 0.1F;
    private int _readstd = 2;
    private int _readdepth = 10;
    Normal normalDist = new Normal(0, 1);
    SNPs _snps;
    int nrSNPs;

    public RawDataSingleSubclone(CopyNumbers copyNumbers, float purity, float ploidy)
    {
        _copynumbers = copyNumbers;
        _purity = purity;
        _ploidy = ploidy;

        // Create SNPs
        int nrSNPsInitial = 100;
        _snps = new SNPs(copyNumbers.IsFemale, nrSNPsInitial);
        nrSNPs = _snps.Length;

        _baf = new float[nrSNPs];
        _logr = new float[nrSNPs];
        CalcRawData();
    }

    private void CalcRawData()
    {
        double randomGaussianValue = normalDist.Sample();

        int curSegmentId = 0;
        var segments = _copynumbers.GetAllSegments();
        var copyNumbers = _copynumbers.GetAllCopyNumbers();
        var snps = _snps.AllSNPs;
        for (int i = 0; i < _snps.Length; i++)
        {
            // Sort snps
            while (snps[i].AbsPos > ReferenceGenome.ChromosomeAbsoluteStart(segments[curSegmentId].ChromId.ChromNum) + segments[curSegmentId].End)
            {
                curSegmentId++;
            }
            if ((copyNumbers[curSegmentId].Item1 + copyNumbers[curSegmentId].Item2) > 0)
            {
                int totalReads = (int)Math.Max(0, Math.Round((copyNumbers[curSegmentId].Item1 + copyNumbers[curSegmentId].Item2) * _readdepth + normalDist.Sample() * _readstd));


                int readsH1 = (int)Math.Max(0, Math.Round((copyNumbers[curSegmentId].Item1) * _readdepth + normalDist.Sample() * _readstd));
                int readsH2 = (int)Math.Max(0, Math.Round((copyNumbers[curSegmentId].Item2) * _readdepth + normalDist.Sample() * _readstd));

                // _logr[i] = (copyNumbers[curSegmentId].Item1 + copyNumbers[curSegmentId].Item2) * _readdepth;
                _logr[i] = (float)Math.Log2((readsH1 + readsH2) / (2f * _readdepth));

                if (snps[i].Heterozygous)
                {
                    _baf[i] = ((readsH1 + readsH2) > 0) ? (float)readsH2 / ((float)readsH1 + (float)readsH2) : -1;
                }
                else
                {
                    _baf[i] = 0f;
                }
                // TODO: Extra noise for logr and baf
            }
        }
    }
    public string BAFToTSV(bool printHeader = false)
    {
        var snps = _snps.AllSNPs;
        var outputString = printHeader ? "\tchrom\tpos\tBAF\n" : "";
        for (int i = 0; i < nrSNPs; i++) {
            outputString += $"snpID\t{snps[i].Chrom}\t{snps[i].Pos}\t{_baf[i]}\n";
        }
        return outputString;
    }
    public string logRToTSV(bool printHeader = false)
    {
        var snps = _snps.AllSNPs;
        var outputString = printHeader ? "\tchrom\tpos\tBAF\n" : "";
        for (int i = 0; i < nrSNPs; i++) {
            outputString += $"snpID\t{snps[i].Chrom}\t{snps[i].Pos}\t{_logr[i]}\n";
        }
        return outputString;
    }
}
public struct SNPs
{
    public IEnumerable<ChromNum> ReferenceChromosomes;
    private List<SNP> _snps = new List<SNP>();
    private float _hetrate = 1F;
    private readonly Random _random = new();
    public SNPs(bool isFemale, int nrsnps)
    {
        ReferenceChromosomes = ReferenceGenome.GetChromosomes(isFemale);

        var totGenomeLength = ReferenceGenome.TotalLength(isFemale);

        foreach (var chrom in ReferenceChromosomes)
        {
            var curChromLength = ReferenceGenome.ChromosomeLengthMap[chrom];
            var nrChromSnps = (int)Math.Floor((double) nrsnps * curChromLength / totGenomeLength);
            for (int i = 0; i < nrChromSnps; i++)
            {
                _snps.Add(new SNP(chrom, _random.Next(1, curChromLength - 1), _random.NextDouble() < _hetrate));
            }
        }

        // Sort all SNPs
        _snps = _snps.OrderBy(snp => snp.AbsPos).ToList();
    }

    public override string ToString()
        => string.Join("\n", _snps.Select(r => r.ToString()));

    public List<SNP> AllSNPs => _snps.ToList();

    public int Length => _snps.Count();

}


public struct SNP
{
    public ChromNum Chrom;
    public int Pos;
    public long AbsPos;
    public Nucleotides Ref;
    public Nucleotides Alt;
    public bool Heterozygous;

    public SNP(ChromNum chrom, int pos, bool heterozygous, Nucleotides re, Nucleotides alt)
    {
        Chrom = chrom;
        Pos = pos;
        AbsPos = ReferenceGenome.ChromosomeAbsoluteStart(chrom) + pos;
        Ref = re;
        Alt = alt;
        Heterozygous = heterozygous;
    }
    public SNP(ChromNum chrom, int pos, bool heterozygous)
    {
        Chrom = chrom;
        Pos = pos;
        AbsPos = ReferenceGenome.ChromosomeAbsoluteStart(chrom) + pos;
        Ref = Nucleotides.A;
        Alt = Nucleotides.G;
        Heterozygous = heterozygous;
    }

    public override string ToString()
        => $"{Chrom} {Pos}: {Ref} / {Alt} ({(Heterozygous ? "het" : "hom")})";
}

public enum Nucleotides
{
    G, T, C, A
}
