// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com
using System.Text;

namespace SimChA.DataTypes;

public class GenRef
{
    public string Name { get; }
    public Dictionary<string, int> ChrLengths { get; }
    public Dictionary<string, SexEnum> ChrSex { get; }
    public int AutosomeCount { get; }
    public int ChrCount { get; }
    private Region[] XYGenome { get; }
    private Region[] XXGenome { get; }
    private Region[] Autosomes { get; }
    
    private long XYLinLen { get; }
    private long XXLinLen { get; }
    private long XYGenomeLen { get; }
    private long XXGenomeLen { get; }
    private List<string> XYChrs { get; }
    private List<string> XXChrs { get; }
    public List<string> AllChrs { get; }
    public List<string> AutosomeChrs { get; }
    
    public string YChrName { get; }
    
    public string XChrName { get; }
    
    public Dictionary<string, StringBuilder>? GenContentsDict { get; set;}

    public Region[] GetGenotype(bool sexXX)
        => IncludeSexChromosomes ? (sexXX ? XXGenome : XYGenome) : Autosomes;
    
    public long GetGenomeLen(bool sexXX, bool diploid = true)
    {
        if (!IncludeSexChromosomes)
        {
            return AutosomeLen;
        }
        else 
        {
            return diploid ? (sexXX ? XXGenomeLen : XYGenomeLen) : (sexXX ? XXLinLen : XYLinLen);
        }
    }
    
    public IEnumerable<string> ChrIDsForSex(bool sexXX)
        => IncludeSexChromosomes ? (sexXX ? XXChrs : XYChrs) : AutosomeChrs;
    
    public long AutosomeLen {get;}
    public bool IncludeSexChromosomes { get; }

    public Dictionary<GeneListType, Dictionary<string, List<Gene>>> GeneLists { get; }

    public GenRef(string name, Dictionary<string, int> chrLengths, Dictionary<string, SexEnum> chrSex, 
        Dictionary<GeneListType, Dictionary<string, List<Gene>>> geneList, bool includeSexChromosomes, Dictionary<string, StringBuilder>? genContentsDict = null)
    {
        Name = name;
        ChrLengths = chrLengths;
        ChrSex = chrSex;
        AutosomeCount = chrSex.Count(x => x.Value == SexEnum.Both);
        IncludeSexChromosomes = includeSexChromosomes;
        ChrCount = AutosomeCount * 2 + (chrSex.Count - AutosomeCount);
        XYChrs = chrSex.Select(pair => pair.Key).ToList();
        XXChrs = chrSex.Where(pair => pair.Value != SexEnum.Male).Select(pair => pair.Key).ToList();
        AllChrs = chrSex.Select(pair => pair.Key).ToList();
        YChrName = chrSex.Where(pair => pair.Value == SexEnum.Male).Select(pair => pair.Key).FirstOrDefault("");
        XChrName = chrSex.Where(pair => pair.Value == SexEnum.Female).Select(pair => pair.Key).FirstOrDefault("");
        bool useSNV = genContentsDict != null;
        var haplotypeOneF = CreateHaplotype(true, true, useSNV);
        var haplotypeTwoF = CreateHaplotype(false, true, useSNV);
        var haplotypeOneM = CreateHaplotype(true, false, useSNV);
        var haplotypeTwoM = CreateHaplotype(false, false, useSNV);
        XYGenome = haplotypeOneM.Concat(haplotypeTwoM).ToArray();
        XXGenome = haplotypeOneF.Concat(haplotypeTwoF).ToArray();
        XYLinLen = XYChrs.Select(c => (long) chrLengths[c]).Sum();
        XXLinLen = XXChrs.Select(c => (long) chrLengths[c]).Sum();
        XYGenomeLen = XYGenome.Sum(r => r.Length);
        XXGenomeLen = XXGenome.Sum(r => r.Length);
        AutosomeChrs = chrSex.Where(pair => pair.Value != SexEnum.Male && pair.Value != SexEnum.Female).Select(pair => pair.Key).ToList();
        AutosomeLen = AutosomeChrs.Select(c => (long) chrLengths[c]).Sum();
        var autohaplotypeOne = CreateAutosomeHaplotype(true, useSNV);
        var autohaplotypeTwo = CreateAutosomeHaplotype(false, useSNV);
        Autosomes = autohaplotypeOne.Concat(autohaplotypeTwo).ToArray();
        GeneLists = geneList;
        GenContentsDict = genContentsDict;
    }

    private Region GetRegion(string chrNo, bool isFirstHaplotype, bool useSNV) => 
        new(0, ChrLengths[chrNo], chrNo, isFirstHaplotype, true, useSNV ? new Dictionary<long, Nucleotide>() : null);

    private IEnumerable<Region> CreateHaplotype(bool isFirstHaplotype, bool isFemale, bool useSNV = false)
    {
        var nonGender = ChrSex.Select(x => x.Key).Where(x => ChrSex[x] == SexEnum.Both);
        var sexChr = ChrSex.Select(x => x.Key).Where(x =>
            isFirstHaplotype | isFemale ? ChrSex[x] == SexEnum.Female : ChrSex[x] == SexEnum.Male);
        var all = nonGender.Concat(sexChr);
        return all.Select(num => GetRegion(num, isFirstHaplotype, useSNV));
    }
    private IEnumerable<Region> CreateAutosomeHaplotype(bool isFirstHaplotype, bool useSNV = false)
    {
        var autosomes = ChrSex.Select(x => x.Key).Where(x => ChrSex[x] == SexEnum.Both);
        return autosomes.Select(num => GetRegion(num, isFirstHaplotype, useSNV));
    }
}