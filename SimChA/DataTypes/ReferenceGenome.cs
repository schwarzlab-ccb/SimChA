// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public static class ReferenceGenome
{
    public static readonly Dictionary<ChromNum, int> ChromosomeLengthMap = new()
    {
        {
            ChromNum.chr1,
            247_249_719
        },
        {
            ChromNum.chr2,
            242_951_149
        },
        {
            ChromNum.chr3,
            199_501_827
        },
        {
            ChromNum.chr4,
            191_273_063
        },
        {
            ChromNum.chr5,
            180_857_866
        },
        {
            ChromNum.chr6,
            170_899_992
        },
        {
            ChromNum.chr7,
            158_821_424
        },
        {
            ChromNum.chr8,
            146_274_826
        },
        {
            ChromNum.chr9,
            140_273_252
        },
        {
            ChromNum.chr10,
            135_374_737
        },
        {
            ChromNum.chr11,
            134_452_384
        },
        {
            ChromNum.chr12,
            132_349_534
        },
        {
            ChromNum.chr13,
            114_142_980
        },
        {
            ChromNum.chr14,
            106_368_585
        },
        {
            ChromNum.chr15,
            100_338_915
        },
        {
            ChromNum.chr16,
            88_827_254
        },
        {
            ChromNum.chr17,
            78_774_742
        },
        {
            ChromNum.chr18,
            76_117_153
        },
        {
            ChromNum.chr19,
            63_811_651
        },
        {
            ChromNum.chr20,
            62_435_964
        },
        {
            ChromNum.chr21,
            46_944_323
        },
        {
            ChromNum.chr22,
            49_691_432
        },
        {
            ChromNum.chrX,
            154_913_754
        },
        {
            ChromNum.chrY,
            57_772_954
        }
    };

    public static readonly Dictionary<ChromNum, long> ChromosomeStartMap;

    private static IEnumerable<Region> CreateHaplotype(bool isFirstHaplotype, bool isFemale)
    {
        var nonGender = Enum.GetValues<ChromNum>().Take(22);
        var sexChr = (isFirstHaplotype | isFemale) ? ChromNum.chrX : ChromNum.chrY;
        var all = nonGender.Concat(new[] { sexChr });
        return all.Select(num => GetRegion(num, isFirstHaplotype));
    }

    static ReferenceGenome()
    {
        var haplotypeOneF = CreateHaplotype(true, true);
        var haplotypeTwoF = CreateHaplotype(false, true);
        GenotypeF = haplotypeOneF.Concat(haplotypeTwoF).ToArray();
        var haplotypeOneM = CreateHaplotype(true, false);
        var haplotypeTwoM = CreateHaplotype(false, false);
        GenotypeM = haplotypeOneM.Concat(haplotypeTwoM).ToArray();

        long start = 0;
        // ChromosomeStartMap = ChromosomeLengthMap.ToDictionary<ChromNum, long>(pair => pair.Key, pair => { start += pair.Value; });
    }

    private static Region[] GenotypeM { get; }
    private static Region[] GenotypeF { get; }

    public static Region[] GetGenotype(bool isFemale) 
        => isFemale ? GenotypeF : GenotypeM;

    public static long TotalLength(bool isFemale)
        => GetChromosomes(isFemale).Select(chrom => (long)ChromosomeLengthMap[chrom]).Sum();

    public static IEnumerable<ChromNum> GetChromosomes(bool isFemale)
        => Enum.GetValues<ChromNum>().Take(22).Append(isFemale ? ChromNum.chrX : ChromNum.chrY);

    public static Region GetRegion(ChromNum chromNum, bool isFirstHaplotype = true) => 
        new(0, ChromosomeLengthMap[chromNum] + 1, new ChromID(chromNum, isFirstHaplotype));

    public static long ChromosomeAbsoluteStart(ChromNum chromNum) 
    {
        var ReferenceChromosomes = Enum.GetValues<ChromNum>();
        long sum = 0;
        foreach (var c in ReferenceChromosomes)
        {
            if (c == chromNum) break;
            sum += ChromosomeLengthMap[c];
        }
        return sum;
    }
}