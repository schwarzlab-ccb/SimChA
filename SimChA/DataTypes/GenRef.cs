// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com
using System.Linq;
using System.Text;

namespace SimChA.DataTypes;

public class GenRef
{
    public string Name { get; }
    public Dictionary<string, int> ChrLengths { get; }
    public Dictionary<string, SexEnum> ChrSex { get; }
    public int AutosomeCount { get; }
    public int ChrCount 
        => IncludeSexChromosomes ? AutosomeCount * 2 + (ChrSex.Count - AutosomeCount) : AutosomeCount * 2;
    private List<List<Region>> XYGenome { get; }
    private List<List<Region>> XXGenome { get; }
    private List<List<Region>> Autosomes { get; }
    private long XYLinLen { get; }
    private long XXLinLen { get; }
    public long AutosomeLinLen { get; }
    private long XYGenomeLen { get; }
    private long XXGenomeLen { get; }
    private List<string> XYChrs { get; }
    private List<string> XXChrs { get; }
    public List<string> AllChrs { get; }
    private List<string> AutosomeChrs { get; }
    
    public string YChrName { get; }
    
    public string XChrName { get; }
    
    public Dictionary<string, StringBuilder>? GenContentsDict { get; set;}
    
    public long GetGenomeLen(bool sexXX, bool diploid = true)
    {
        if (!IncludeSexChromosomes)
        {
            return diploid ? AutosomeLen : AutosomeLinLen;
        }
        else 
        {
            return diploid ? (sexXX ? XXGenomeLen : XYGenomeLen) : (sexXX ? XXLinLen : XYLinLen);
        }
    }
    
    public IEnumerable<string> ChrIDsForSex(bool sexXX)
        => sexXX ? XXChrs : XYChrs;
    
    public IEnumerable<string> ChrIDsForAutosomes()
        => AutosomeChrs;
    
    public long AutosomeLen {get;}
    public bool IncludeSexChromosomes { get; set;}
    

    public List<List<Region>> GetGenotype(bool sexXX)
        => IncludeSexChromosomes ? (sexXX ? XXGenome : XYGenome) : Autosomes;

    public Dictionary<GeneListType, Dictionary<string, List<Gene>>> GeneLists { get; }

    public GenRef(string name, Dictionary<string, int> chrLengths, Dictionary<string, SexEnum> chrSex,
        (Dictionary<string, (int, int)> p, Dictionary<string, (int, int)> q) centromeres,
        Dictionary<GeneListType, Dictionary<string, List<Gene>>> geneList, bool includeSexChromosomes, Dictionary<string, StringBuilder>? genContentsDict = null)
    {
        Name = name;
        ChrLengths = chrLengths;
        ChrSex = chrSex;
        AutosomeCount = chrSex.Count(x => x.Value == SexEnum.Both);
        IncludeSexChromosomes = includeSexChromosomes;
        XYChrs = chrSex.Select(pair => pair.Key).ToList();
        XXChrs = chrSex.Where(pair => pair.Value != SexEnum.Male).Select(pair => pair.Key).ToList();
        AllChrs = chrSex.Select(pair => pair.Key).ToList();
        YChrName = chrSex.Where(pair => pair.Value == SexEnum.Male).Select(pair => pair.Key).FirstOrDefault("");
        XChrName = chrSex.Where(pair => pair.Value == SexEnum.Female).Select(pair => pair.Key).FirstOrDefault("");
        bool useSNV = genContentsDict != null;
        // Create the haplotypes
        var haplotypeOneF = CreateHaplotype(true, true, useSNV);
        var haplotypeTwoF = CreateHaplotype(false, true, useSNV);
        var haplotypeOneM = CreateHaplotype(true, false, useSNV);
        var haplotypeTwoM = CreateHaplotype(false, false, useSNV);
        XYGenome = haplotypeOneM.Concat(haplotypeTwoM).ToList();
        XXGenome = haplotypeOneF.Concat(haplotypeTwoF).ToList();

        XYLinLen = XYChrs.Select(c => (long) chrLengths[c]).Sum();
        XXLinLen = XXChrs.Select(c => (long) chrLengths[c]).Sum();
        XYGenomeLen = XYGenome.Sum(regions => regions.Sum(r => r.Length));
        XXGenomeLen = XXGenome.Sum(regions => regions.Sum(r => r.Length));
        AutosomeChrs = chrSex.Where(pair => pair.Value != SexEnum.Male && pair.Value != SexEnum.Female).Select(pair => pair.Key).ToList();
        AutosomeLinLen = AutosomeChrs.Select(c => (long) chrLengths[c]).Sum();
        var autohaplotypeOne = CreateHaplotype(true, true, true, useSNV);
        var autohaplotypeTwo = CreateHaplotype(true, true, true, useSNV);
        Autosomes = autohaplotypeOne.Concat(autohaplotypeTwo).ToList();
        AutosomeLen = Autosomes.Sum(regions => regions.Sum(r => r.Length));
        GeneLists = geneList;
        GenContentsDict = genContentsDict;
    }

    private Region GetRegion(string chrNo, bool isFirstHaplotype, bool useSNV) => 
        new(0, ChrLengths[chrNo], chrNo, isFirstHaplotype, true, useSNV ? new Dictionary<long, Nucleotide>() : null);

    private IEnumerable<List<Region>> CreateAutosomeHaplotype(bool isFirstHaplotype, bool useSNV = false)
    {
        var autosomes = ChrSex.Select(x => x.Key).Where(x => ChrSex[x] == SexEnum.Both);
        return autosomes.Select(num => new List<Region>{GetRegion(num, isFirstHaplotype, useSNV)});
    }

    private IEnumerable<List<Region>> CreateHaplotype(bool isFirstHaplotype, bool isFemale, bool autosomesOnly = false, bool useSNV = false)
    {
        var chrs = ChrSex.Select(x => x.Key).Where(x => ChrSex[x] == SexEnum.Both);
        if (!autosomesOnly)
        {
            var sexChr = ChrSex.Select(x => x.Key).Where(x =>
                isFirstHaplotype | isFemale ? ChrSex[x] == SexEnum.Female : ChrSex[x] == SexEnum.Male);
            chrs = chrs.Concat(sexChr);
        }
        return chrs.Select(num => new List<Region>{GetRegion(num, isFirstHaplotype, useSNV)});
    }
    
}