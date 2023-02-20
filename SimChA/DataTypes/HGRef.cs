// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com
// Edited by Felix Schifferdecker to match the GRCh38.p13

namespace SimChA.DataTypes;

public static class HGRef
{
    public const int CHR_COUNT = 46;
    public const int AUTOSOME_COUNT = 22;

    private static readonly Dictionary<ChrNo, int> ChromosomeLengthMap = new()
    {
        {
            ChrNo.chr1,
            249_956_942
        },
        {
            ChrNo.chr2,
            242_508_799
        },
        {
            ChrNo.chr3,
            198_450_956
        },
        {
            ChrNo.chr4,
            190_424_264
        },
        {
            ChrNo.chr5,
            181_630_948
        },
        {
            ChrNo.chr6,
            170_805_979
        },
        {
            ChrNo.chr7,
            159_345_973
        },
        {
            ChrNo.chr8,
            145_138_636
        },
        {
            ChrNo.chr9,
            138_688_728
        },
        {
            ChrNo.chr10,
            133_797_422
        },
        {
            ChrNo.chr11,
            135_186_938
        },
        {
            ChrNo.chr12,
            133_275_309
        },
        {
            ChrNo.chr13,
            114_364_328
        },
        {
            ChrNo.chr14,
            108_136_338
        },
        {
            ChrNo.chr15,
            102_439_437
        },
        {
            ChrNo.chr16,
            92_211_104
        },
        {
            ChrNo.chr17,
            83_836_422
        },
        {
            ChrNo.chr18,
            80_373_285
        },
        {
            ChrNo.chr19,
            58_617_616
        },
        {
            ChrNo.chr20,
            64_444_167
        },
        {
            ChrNo.chr21,
            46_709_983
        },
        {
            ChrNo.chr22,
            51_857_516
        },
        {
            ChrNo.chrX,
            156_040_895
        },
        {
            ChrNo.chrY,
            57_264_655
        }
    };
    
    private static long maleGenomeLen;
    private static long femaleGenomeLen;

    static HGRef()
    {
        var haplotypeOneF = CreateHaplotype(true, true);
        var haplotypeTwoF = CreateHaplotype(false, true);
        GenotypeF = haplotypeOneF.Concat(haplotypeTwoF).ToArray();
        var haplotypeOneM = CreateHaplotype(true, false);
        var haplotypeTwoM = CreateHaplotype(false, false);
        GenotypeM = haplotypeOneM.Concat(haplotypeTwoM).ToArray();

        femaleGenomeLen = TotalLength(true);
        maleGenomeLen = TotalLength(false);
    }

    private static Region[] GenotypeM { get; }
    private static Region[] GenotypeF { get; }

    private static IEnumerable<Region> CreateHaplotype(bool isFirstHaplotype, bool isFemale)
    {
        var nonGender = Enum.GetValues<ChrNo>().Take(22);
        var sexChr = isFirstHaplotype | isFemale ? ChrNo.chrX : ChrNo.chrY;
        var all = nonGender.Concat(new[] { sexChr });
        return all.Select(num => GetRegion(num, isFirstHaplotype));
    }

    public static Region[] GetGenotype(bool isFemale)
        => isFemale ? GenotypeF : GenotypeM;

    public static long TotalLength(bool isFemale)
    {
        long length = 2 * Enum.GetValues<ChrNo>().Take(22).Select(chrom => (long)ChromosomeLengthMap[chrom]).Sum();
        if (isFemale)
        {
            length += 2 * ChromosomeLengthMap[ChrNo.chrX];
        }
        else
        {
            length += ChromosomeLengthMap[ChrNo.chrX] + ChromosomeLengthMap[ChrNo.chrY];
        }
        return length;
    }

    public static IEnumerable<ChrNo> ChrIDsForSex(bool isFemale)
        => Enum.GetValues<ChrNo>().Take(isFemale ? AUTOSOME_COUNT + 1: AUTOSOME_COUNT + 2);

    public static Region GetRegion(ChrNo chrNo, bool isFirstHaplotype = true) =>
        new(0, ChromosomeLengthMap[chrNo], new ChrID(chrNo, isFirstHaplotype));
    
    public static long GetGenomeSize(bool isFemale)
        => isFemale ? femaleGenomeLen : maleGenomeLen;
}