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

    private static IEnumerable<Region> CreateHaplotype(bool isHaplotypeOne)
        => Enum.GetValues<ChromNum>()
            .Select(num
                => new Region(0, ChromosomeLengthMap[num] + 1, new ChromID(num, isHaplotypeOne))
            );

    static ReferenceGenome()
    {
        var haplotypeOne = CreateHaplotype(true);
        var haplotypeTwo = CreateHaplotype(false);
        Genotype = haplotypeOne.Concat(haplotypeTwo).ToArray();
    }

    public static Region[] Genotype { get; }
}