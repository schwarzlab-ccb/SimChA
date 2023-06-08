// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public static class HGRef
{
    public const int CHR_COUNT = 46;
    public const int AUTOSOME_COUNT = 22;
    
    public static string Sex(bool isXX) => isXX ? "XX" : "XY";

    // https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37
    private static readonly Dictionary<ChrNo, int> HG19Chrs = new()
    {
        {ChrNo.chr1, 249_250_621},
        {ChrNo.chr2, 243_199_373},
        {ChrNo.chr3, 198_022_430},
        {ChrNo.chr4, 191_154_276},
        {ChrNo.chr5, 180_915_260},
        {ChrNo.chr6, 171_115_067},
        {ChrNo.chr7, 159_138_663},
        {ChrNo.chr8, 146_364_022},
        {ChrNo.chr9, 141_213_431},
        {ChrNo.chr10, 135_534_747},
        {ChrNo.chr11, 135_006_516},
        {ChrNo.chr12, 133_851_895},
        {ChrNo.chr13, 115_169_878},
        {ChrNo.chr14, 107_349_540},
        {ChrNo.chr15, 102_531_392},
        {ChrNo.chr16, 90_354_753},
        {ChrNo.chr17, 81_195_210},
        {ChrNo.chr18, 78_077_248},
        {ChrNo.chr19, 59_128_983},
        {ChrNo.chr20, 63_025_520},
        {ChrNo.chr21, 48_129_895},
        {ChrNo.chr22, 51_304_566},
        {ChrNo.chrX, 155_270_560},
        {ChrNo.chrY, 59_373_566}
    };

    // https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38
    private static readonly Dictionary<ChrNo, int> HG38Chrs = new()
    {
        {ChrNo.chr1, 248_956_422},
        {ChrNo.chr2, 242_193_529},
        {ChrNo.chr3, 198_295_559},
        {ChrNo.chr4, 190_214_555},
        {ChrNo.chr5, 181_538_259}, 
        {ChrNo.chr6, 170_805_979},
        {ChrNo.chr7, 159_345_973},
        {ChrNo.chr8, 145_138_636},
        {ChrNo.chr9, 138_394_717},
        {ChrNo.chr10, 133_797_422},
        {ChrNo.chr11, 135_086_622},
        {ChrNo.chr12, 133_275_309},
        {ChrNo.chr13, 114_364_328},
        {ChrNo.chr14, 107_043_718},
        {ChrNo.chr15, 101_991_189},
        {ChrNo.chr16, 90_338_345},
        {ChrNo.chr17, 83_257_441},
        {ChrNo.chr18, 80_373_285},
        {ChrNo.chr19, 58_617_616},
        {ChrNo.chr20, 64_444_167},
        {ChrNo.chr21, 46_709_983},
        {ChrNo.chr22, 50_818_468},
        {ChrNo.chrX, 156_040_895},
        {ChrNo.chrY, 57_227_415}
    };

    private static GenomeAssembly _assembly;

    private static Dictionary<ChrNo, int> GetChrs(GenomeAssembly assembly) 
        => assembly == GenomeAssembly.hg19 ? HG19Chrs : HG38Chrs;

    public static GenomeAssembly Assembly
    {
        get
        {
            if (_assembly == GenomeAssembly.none)
            {
                throw new Exception("Genome assembly not set.");
            }
            return _assembly;
        }
        set => _assembly = value;
    }

    private static Dictionary<GenomeAssembly, GenomeReference> Refs { get; }

    public static Region[] GetGenotype(bool sexXX)
        => sexXX ? Refs[Assembly].XXGenome : Refs[Assembly].XYGenome;

    public static long GetGenomeLen(bool sexXX)
        => sexXX ? Refs[Assembly].XXGenomeLen : Refs[Assembly].XYGenomeLen;

    public static long GetChromLen(ChrNo chrNo)
        => GetChrs(Assembly)[chrNo];
    
    public static IEnumerable<ChrNo> ChrIDsForSex(bool sexXX)
        => Enum.GetValues<ChrNo>().Take(sexXX ? AUTOSOME_COUNT + 1 : AUTOSOME_COUNT + 2);
    
    public static Region GetRegion(ChrNo chrNo, bool isFirstHaplotype = true) 
        => new(0, GetChrs(Assembly)[chrNo], chrNo, isFirstHaplotype);

    static HGRef()
    {
        Refs = new Dictionary<GenomeAssembly, GenomeReference>();
        foreach (var assembly in Enum.GetValues<GenomeAssembly>())
        {
            var haplotypeOneF = CreateHaplotype(true, true, assembly);
            var haplotypeTwoF = CreateHaplotype(false, true, assembly);
            var haplotypeOneM = CreateHaplotype(true, false, assembly);
            var haplotypeTwoM = CreateHaplotype(false, false, assembly);
            var xyGenome = haplotypeOneM.Concat(haplotypeTwoM).ToArray();
            var xxGenome = haplotypeOneF.Concat(haplotypeTwoF).ToArray();
            long xyGenomeLen = xyGenome.Sum(r => r.Length);
            long xxGenomeLen = xxGenome.Sum(r => r.Length);
            Refs[assembly] = new GenomeReference(xyGenome, xyGenomeLen, xxGenome, xxGenomeLen);
        }
    }
    
    private static Region GetRegion(ChrNo chrNo, GenomeAssembly genomeAssembly, bool isFirstHaplotype = true) 
        => new(0, GetChrs(genomeAssembly)[chrNo], chrNo, isFirstHaplotype);

    private static IEnumerable<Region> CreateHaplotype(bool isFirstHaplotype, bool isFemale, GenomeAssembly assembly)
    {
        var nonGender = Enum.GetValues<ChrNo>().Take(22);
        var sexChr = isFirstHaplotype | isFemale ? ChrNo.chrX : ChrNo.chrY;
        var all = nonGender.Concat(new[] {sexChr});
        return all.Select(num => GetRegion(num, assembly, isFirstHaplotype));
    }
}