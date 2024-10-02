// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using System.Collections.Immutable;
using System.Text;

namespace SimChA.DataTypes;

public class GenRef
{
    public GenRef(
        string name,
        Dictionary<string, int> chrLengths,
        Dictionary<string, SexEnum> chrSex,
        ImmutableDictionary<string, (long start, long end)> centromeres,
        Dictionary<GeneListType, Dictionary<string, List<Gene>>> geneList,
        bool includeSexChromosomes,
        Dictionary<string, StringBuilder>? genContentsDict = null)
    {
        Name = name;
        ChrLengths = chrLengths;
        ChrSex = chrSex;
        Centromeres = centromeres;
        AutosomesCount = chrSex.Count(x => x.Value == SexEnum.None);
        IncludeSexChromosomes = includeSexChromosomes;
        YChrs = chrSex.Where(pair => pair.Value != SexEnum.Female).Select(pair => pair.Key).ToList();
        XChrs = chrSex.Where(pair => pair.Value != SexEnum.Male).Select(pair => pair.Key).ToList();
        AutChrs = chrSex.Where(pair => pair.Value == SexEnum.None).Select(pair => pair.Key).ToList();
        AllChrs = chrSex.Select(pair => pair.Key).ToList();
        YChrName = chrSex.Where(pair => pair.Value == SexEnum.Male).Select(pair => pair.Key).FirstOrDefault("");
        XChrName = chrSex.Where(pair => pair.Value == SexEnum.Female).Select(pair => pair.Key).FirstOrDefault("");
        bool useSNV = genContentsDict != null;

        // Create the haplotypes
        var haplotypeOneF = CreateHaplotype(SexEnum.Female, true, useSNV);
        var haplotypeTwoF = CreateHaplotype(SexEnum.Female,false, useSNV);
        var haplotypeOneM = CreateHaplotype(SexEnum.Male,true, useSNV);
        var haplotypeTwoM = CreateHaplotype(SexEnum.Male,false, useSNV);
        var haplotypeOneA = CreateHaplotype(SexEnum.None,true, useSNV);
        var haplotypeTwoA = CreateHaplotype(SexEnum.None,false, useSNV);
        XYGenome = haplotypeOneM.Concat(haplotypeTwoM).ToList();
        XXGenome = haplotypeOneF.Concat(haplotypeTwoF).ToList();
        Autosome = haplotypeOneA.Concat(haplotypeTwoA).ToList();

        XYLinLen = YChrs.Select(c => (long)chrLengths[c]).Sum();
        XXLinLen = XChrs.Select(c => (long)chrLengths[c]).Sum();
        AutosomeLinLen = AutChrs.Select(c => (long)chrLengths[c]).Sum();
        XYGenomeLen = XYGenome.Sum(regions => regions.Sum(r => r.Length));
        XXGenomeLen = XXGenome.Sum(regions => regions.Sum(r => r.Length));
        AutosomeLen = Autosome.Sum(regions => regions.Sum(r => r.Length));

        GeneLists = geneList;
        GenContentsDict = genContentsDict;
    }
    public string Name { get; }
    public Dictionary<string, int> ChrLengths { get; }
    public Dictionary<string, SexEnum> ChrSex { get; }
    public ImmutableDictionary<string, (long start, long end)> Centromeres { get; }
    public int AutosomesCount { get; }
    
    public int ChrCount(SexEnum sex, bool diploid = true)
        => (diploid, sex) switch
        {
            (true, SexEnum.None) => AutosomesCount * 2,
            (true, SexEnum.Female) => 46,
            (true, SexEnum.Male) => 46,
            (false, SexEnum.None) => AutosomesCount,
            (false, SexEnum.Female) => 23,
            (false, SexEnum.Male) => 23,
            _ => throw new ArgumentOutOfRangeException($"Missing chromosome counts for {sex}, {diploid}")
        };

    private List<List<Region>> XYGenome { get; }
    private List<List<Region>> XXGenome { get; }
    private List<List<Region>> Autosome { get; }
    private long XYLinLen { get; }
    private long XXLinLen { get; }
    public long AutosomeLinLen { get; }
    private long XYGenomeLen { get; }
    private long XXGenomeLen { get; }
    public long AutosomeLen { get; }
    private List<string> YChrs { get; }
    private List<string> XChrs { get; }
    public List<string> AllChrs { get; }
    private List<string> AutChrs { get; }
    public string YChrName { get; }
    public string XChrName { get; }
    public bool IncludeSexChromosomes { get; set; } // TODO: This should not be part of a reference

    public Dictionary<string, StringBuilder>? GenContentsDict { get; set; }

    public Dictionary<GeneListType, Dictionary<string, List<Gene>>> GeneLists { get; }

    public long GetGenomeLen(SexEnum sex, bool diploid = true)
        => (diploid, sex) switch
        {
            (true, SexEnum.Female) => XXGenomeLen,
            (true, SexEnum.Male) => XYGenomeLen,
            (true, SexEnum.None) => AutosomeLen,
            (false, SexEnum.Female) => XXLinLen,
            (false, SexEnum.Male) => XYLinLen,
            (false, SexEnum.None) => AutosomeLinLen,
            _ => throw new ArgumentOutOfRangeException($"Missing genome length for {sex}, {diploid}")
        };
    
    public IEnumerable<string> ChrIDsForSex(SexEnum sex)
        => sex switch
        {
            SexEnum.Female => XChrs, 
            SexEnum.Male => AllChrs, 
            _ => AutChrs
        };
    
    public List<string> ChrIDsForHap(SexEnum sex, bool diploid = true)
        => (diploid, sex) switch
        {
            (true, SexEnum.Female) => XChrs,
            (true, SexEnum.Male) => XChrs,
            (true, SexEnum.None) => AutChrs,
            (false, SexEnum.Female) => XChrs,
            (false, SexEnum.Male) => YChrs,
            (false, SexEnum.None) => AutChrs,
            _ => throw new ArgumentOutOfRangeException($"Missing chr IDs for {sex}, {diploid}")
        };

    
    public IEnumerable<string> ChrIDsForAutosomes()
        => ChrIDsForSex(SexEnum.None);

    public IEnumerable<List<Region>> GetGenotype(SexEnum sexEnum)
        => sexEnum switch
        {
            SexEnum.Female => XXGenome,
            SexEnum.Male => XYGenome,
            _ => Autosome
        };
    
    private Region GetRegion(string chrNo, bool isFirstHaplotype, bool useSNV) 
        => new(0, ChrLengths[chrNo], chrNo, isFirstHaplotype, true, useSNV ? new Dictionary<long, Nucleotide>() : null);
    
    private IEnumerable<List<Region>> CreateHaplotype(SexEnum sex, bool firstHap, bool useSNV = false)
        => ChrIDsForHap(sex, firstHap).Select(chr => new List<Region> { GetRegion(chr, firstHap, useSNV) });
}