using System.Collections.Immutable;
using System.Text;

namespace SimChA.Data;

public class GenRef
{
    public GenRef(
        string name,
        Dictionary<string, int> chrLengths,
        Dictionary<string, SexType> chrSex,
        ImmutableDictionary<string, (long start, long end)> centromeres,
        Dictionary<GeneListType, Dictionary<string, List<Gene>>> geneList,
        Dictionary<string, StringBuilder>? genContentsDict = null)
    {
        Name = name;
        ChrLengths = chrLengths;
        ChrSex = chrSex;
        Centromeres = centromeres;
        AutosomesCount = chrSex.Count(x => x.Value == SexType.Any);
        YChrs = chrSex.Where(pair => pair.Value != SexType.Female).Select(pair => pair.Key).ToList();
        XChrs = chrSex.Where(pair => pair.Value != SexType.Male).Select(pair => pair.Key).ToList();
        AutChrs = chrSex.Where(pair => pair.Value == SexType.Any).Select(pair => pair.Key).ToList();
        AllChrs = chrSex.Select(pair => pair.Key).ToList();
        YChrName = chrSex.Where(pair => pair.Value == SexType.Male).Select(pair => pair.Key).FirstOrDefault("");
        XChrName = chrSex.Where(pair => pair.Value == SexType.Female).Select(pair => pair.Key).FirstOrDefault("");
        bool useSNV = genContentsDict != null;

        // Create the haplotypes
        var haplotypeOneF = CreateHaplotype(SexType.Female, true, useSNV);
        var haplotypeTwoF = CreateHaplotype(SexType.Female,false, useSNV);
        var haplotypeOneM = CreateHaplotype(SexType.Male,true, useSNV);
        var haplotypeTwoM = CreateHaplotype(SexType.Male,false, useSNV);
        var haplotypeOneA = CreateHaplotype(SexType.Any,true, useSNV);
        var haplotypeTwoA = CreateHaplotype(SexType.Any,false, useSNV);
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
    public Dictionary<string, SexType> ChrSex { get; }
    public ImmutableDictionary<string, (long start, long end)> Centromeres { get; }
    public int AutosomesCount { get; }
    
    public int ChrCount(SexType sex, bool diploid = true)
        => (diploid, sex) switch
        {
            (true, SexType.Female) => XChrs.Count*2,
            (true, SexType.Male) => YChrs.Count*2,
            (true, SexType.Any) => AutChrs.Count*2,
            (false, SexType.Female) => XChrs.Count,
            (false, SexType.Male) => AllChrs.Count,
            (false, SexType.Any) => AutChrs.Count,
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

    public Dictionary<string, StringBuilder>? GenContentsDict { get; set; }

    public Dictionary<GeneListType, Dictionary<string, List<Gene>>> GeneLists { get; }

    public long GetGenomeLen(SexType sex, bool diploid = true)
        => (diploid, sex) switch
        {
            (true, SexType.Female) => XXGenomeLen,
            (true, SexType.Male) => XYGenomeLen,
            (true, SexType.Any) => AutosomeLen,
            (false, SexType.Female) => XXLinLen,
            (false, SexType.Male) => XYLinLen,
            (false, SexType.Any) => AutosomeLinLen,
            _ => throw new ArgumentOutOfRangeException($"Missing genome length for {sex}, {diploid}")
        };
    
    public IEnumerable<string> ChrIDsForSex(SexType sex)
        => sex switch
        {
            SexType.Female => XChrs, 
            SexType.Male => AllChrs, 
            _ => AutChrs
        };
    
    public List<string> ChrIDsForHap(SexType sex, bool firstHaplotype = true)
        => (firstHaplotype, sex) switch
        {
            (true, SexType.Female) => XChrs,
            (true, SexType.Male) => XChrs,
            (true, SexType.Any) => AutChrs,
            (false, SexType.Female) => XChrs,
            (false, SexType.Male) => YChrs,
            (false, SexType.Any) => AutChrs,
            _ => throw new ArgumentOutOfRangeException($"Missing chr IDs for {sex}, {firstHaplotype}")
        };

    
    public IEnumerable<string> ChrIDsForAutosomes()
        => ChrIDsForSex(SexType.Any);

    public IEnumerable<List<Region>> GetGenotype(SexType sexType)
        => sexType switch
        {
            SexType.Female => XXGenome,
            SexType.Male => XYGenome,
            _ => Autosome
        };
    
    private Region GetRegion(string chrNo, bool isFirstHaplotype, bool useSNV) 
        => new(0, ChrLengths[chrNo], chrNo, isFirstHaplotype, useSNV ? new List<SNV>() : null);
    
    private IEnumerable<List<Region>> CreateHaplotype(SexType sex, bool firstHap, bool useSNV = false)
        => ChrIDsForHap(sex, firstHap).Select(chr => new List<Region> { GetRegion(chr, firstHap, useSNV) });
}