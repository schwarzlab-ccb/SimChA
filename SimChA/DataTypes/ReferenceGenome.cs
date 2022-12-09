// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com
// Edited by Felix Schifferdecker to match the GRCh38.p13

namespace SimChA.DataTypes;

public static class ReferenceGenome
{
    public static readonly Dictionary<ChromNum, int> ChromosomeLengthMap = new()
    {
        {
            ChromNum.chr1,
            249_956_942
        },
        {
            ChromNum.chr2,
            242_508_799
        },
        {
            ChromNum.chr3,
            198_450_956
        },
        {
            ChromNum.chr4,
            190_424_264
        },
        {
            ChromNum.chr5,
            181_630_948
        },
        {
            ChromNum.chr6,
            170_805_979
        },
        {
            ChromNum.chr7,
            159_345_973
        },
        {
            ChromNum.chr8,
            145_138_636
        },
        {
            ChromNum.chr9,
            138_688_728
        },
        {
            ChromNum.chr10,
            133_797_422
        },
        {
            ChromNum.chr11,
            135_186_938
        },
        {
            ChromNum.chr12,
            133_275_309
        },
        {
            ChromNum.chr13,
            114_364_328
        },
        {
            ChromNum.chr14,
            108_136_338
        },
        {
            ChromNum.chr15,
            102_439_437
        },
        {
            ChromNum.chr16,
            92_211_104
        },
        {
            ChromNum.chr17,
            83_836_422
        },
        {
            ChromNum.chr18,
            80_373_285
        },
        {
            ChromNum.chr19,
            58_617_616
        },
        {
            ChromNum.chr20,
            64_444_167
        },
        {
            ChromNum.chr21,
            46_709_983
        },
        {
            ChromNum.chr22,
            51_857_516
        },
        {
            ChromNum.chrX,
            156_040_895
        },
        {
            ChromNum.chrY,
            57_264_655
        }
    };

    public static readonly Dictionary<ChromNum, long> ChromosomeStartMap;

    static ReferenceGenome()
    {
        var haplotypeOneF = CreateHaplotype(true, true);
        var haplotypeTwoF = CreateHaplotype(false, true);
        GenotypeF = haplotypeOneF.Concat(haplotypeTwoF).ToArray();
        var haplotypeOneM = CreateHaplotype(true, false);
        var haplotypeTwoM = CreateHaplotype(false, false);
        GenotypeM = haplotypeOneM.Concat(haplotypeTwoM).ToArray();

        long start = -ChromosomeLengthMap.First().Value; // Has to be subtracted, otherwise ends are computed
        ChromosomeStartMap = ChromosomeLengthMap.ToDictionary(pair => pair.Key, pair => start += pair.Value);
    }

    private static Region[] GenotypeM { get; }
    private static Region[] GenotypeF { get; }

    private static IEnumerable<Region> CreateHaplotype(bool isFirstHaplotype, bool isFemale)
    {
        var nonGender = Enum.GetValues<ChromNum>().Take(22);
        var sexChr = isFirstHaplotype | isFemale ? ChromNum.chrX : ChromNum.chrY;
        var all = nonGender.Concat(new[] { sexChr });
        return all.Select(num => GetRegion(num, isFirstHaplotype));
    }

    public static Region[] GetGenotype(bool isFemale)
        => isFemale ? GenotypeF : GenotypeM;

    public static long TotalLength(bool isFemale)
        => GetChromosomes(isFemale).Select(chrom => (long)ChromosomeLengthMap[chrom]).Sum();

    public static IEnumerable<ChromNum> GetChromosomes(bool isFemale)
        => Enum.GetValues<ChromNum>().Take(22).Append(isFemale ? ChromNum.chrX : ChromNum.chrY);

    public static Region GetRegion(ChromNum chromNum, bool isFirstHaplotype = true) =>
        new(0, ChromosomeLengthMap[chromNum] + 1, new ChromID(chromNum, isFirstHaplotype));
}