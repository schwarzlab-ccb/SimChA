// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

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
    private long XYGenomeLen { get; }
    private long XXGenomeLen { get; }
    private List<string> XYChrs { get; }
    private List<string> XXChrs { get; }
    public List<string> AllChrs { get; }
    
    public string YChrName { get; }
    
    public string XChrName { get; }
    
    public List<GenContents>? GenContentsList { get; }

    public Region[] GetGenotype(bool sexXX)
        => sexXX ? XXGenome : XYGenome;
    
    public long GetGenomeLen(bool sexXX)
        => sexXX ? XXGenomeLen : XYGenomeLen;
    
    public IEnumerable<string> ChrIDsForSex(bool sexXX)
        => sexXX ? XXChrs : XYChrs;

    public Dictionary<GeneListType, Dictionary<string, List<Gene>>> GeneLists { get; }

    public GenRef(string name, Dictionary<string, int> chrLengths, Dictionary<string, SexEnum> chrSex, 
        Dictionary<GeneListType, Dictionary<string, List<Gene>>> geneList, List<GenContents>? genContentsList = null)
    {
        Name = name;
        ChrLengths = chrLengths;
        ChrSex = chrSex;
        AutosomeCount = chrSex.Count(x => x.Value == SexEnum.Both);
        ChrCount = AutosomeCount * 2 + (chrSex.Count - AutosomeCount);
        XYChrs = chrSex.Select(pair => pair.Key).ToList();
        XXChrs = chrSex.Where(pair => pair.Value != SexEnum.Male).Select(pair => pair.Key).ToList();
        AllChrs = chrSex.Select(pair => pair.Key).ToList();
        YChrName = chrSex.Where(pair => pair.Value == SexEnum.Male).Select(pair => pair.Key).FirstOrDefault("");
        XChrName = chrSex.Where(pair => pair.Value == SexEnum.Female).Select(pair => pair.Key).FirstOrDefault("");
        bool useSNV = genContentsList != null;
        var haplotypeOneF = CreateHaplotype(true, true, useSNV);
        var haplotypeTwoF = CreateHaplotype(false, true, useSNV);
        var haplotypeOneM = CreateHaplotype(true, false, useSNV);
        var haplotypeTwoM = CreateHaplotype(false, false, useSNV);
        XYGenome = haplotypeOneM.Concat(haplotypeTwoM).ToArray();
        XXGenome = haplotypeOneF.Concat(haplotypeTwoF).ToArray();
        XYGenomeLen = XYGenome.Sum(r => r.Length);
        XXGenomeLen = XXGenome.Sum(r => r.Length);
        GeneLists = geneList;
        GenContentsList = genContentsList;
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
}