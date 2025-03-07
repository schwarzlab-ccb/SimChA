using System.Text;

namespace SimChA.Data;

public class GenRef
{
    public string Name { get; }
    public Dictionary<string, int> ChrLengths { get; }
    public Dictionary<string, SexType> ChrSex { get; }
    public Dictionary<string, GenRange> Centromeres { get; }
    public int AutosomesCount { get; }
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
    public List<string> AutChrs { get; }
    public string YChrName { get; }
    public string XChrName { get; }
    
    private Dictionary<string, StringBuilder>? GenContentsDict { get; }
    public char[] GetGenContents(string chrom, int start, int length)
        => GenContentsDict != null 
            ? GenContentsDict[chrom].ToString(start, length).ToCharArray() 
            : ['N'];
    static public Nucleotide CharToNucleotide(char c) 
    => char.ToUpper(c) switch
        {
            'A' => Nucleotide.A,
            'C' => Nucleotide.C,
            'G' => Nucleotide.G,
            'T' => Nucleotide.T,
            _ => Nucleotide.N  // Default case for invalid characters
        };
    public Nucleotide GetRefBase(string chrom, int location)
        => GenContentsDict != null 
            ? CharToNucleotide(GenContentsDict[chrom].ToString(location, 1)[0])
            : Nucleotide.N;
    public Nucleotide GetRefBaseFromSeq(IEnumerable<string> seq, int location)
    {
        if (GenContentsDict == null) {
            return Nucleotide.N;
        }
        return CharToNucleotide(seq.ElementAt((int)location)[0]);
    }

    public Dictionary<GeneListType, Dictionary<string, List<Gene>>> GeneLists { get; }

    public IEnumerable<Gene> GetTSGs
        => GeneLists[GeneListType.TumorSuppressor].SelectMany(kvp => kvp.Value);
    public IEnumerable<Gene> GetOGs
        => GeneLists[GeneListType.Oncogene].SelectMany(kvp => kvp.Value);
    public IEnumerable<Gene> GetEssGenes
        => GeneLists[GeneListType.Essentiality].SelectMany(kvp => kvp.Value);
    
    public GenRef(
        string name,
        Dictionary<string, int> chrLengths,
        Dictionary<string, SexType> chrSex,
        Dictionary<string, GenRange>  centromeres,
        Dictionary<GeneListType, Dictionary<string, List<Gene>>> geneList,
        Dictionary<string, StringBuilder>? genContentsDict = null)
    {
        Name = name;
        ChrLengths = chrLengths;
        ChrSex = chrSex;
        Centromeres = centromeres;
        GeneLists = geneList;

        AutosomesCount = ChrSex.Count(x => x.Value == SexType.Any);
        YChrs = ChrSex.Where(pair => pair.Value != SexType.Female).Select(pair => pair.Key).ToList();
        XChrs = ChrSex.Where(pair => pair.Value != SexType.Male).Select(pair => pair.Key).ToList();
        AutChrs = ChrSex.Where(pair => pair.Value == SexType.Any).Select(pair => pair.Key).ToList();
        AllChrs = ChrSex.Select(pair => pair.Key).ToList();
        YChrName = ChrSex.Where(pair => pair.Value == SexType.Male).Select(pair => pair.Key).FirstOrDefault("");
        XChrName = ChrSex.Where(pair => pair.Value == SexType.Female).Select(pair => pair.Key).FirstOrDefault("");

        // Create the haplotypes
        var haplotypeOneF = CreateHaplotype(SexType.Female, true);
        var haplotypeTwoF = CreateHaplotype(SexType.Female,false);
        var haplotypeOneM = CreateHaplotype(SexType.Male,true);
        var haplotypeTwoM = CreateHaplotype(SexType.Male,false);
        var haplotypeOneA = CreateHaplotype(SexType.Any,true);
        var haplotypeTwoA = CreateHaplotype(SexType.Any,false);
        XYGenome = haplotypeOneM.Concat(haplotypeTwoM).ToList();
        XXGenome = haplotypeOneF.Concat(haplotypeTwoF).ToList();
        Autosome = haplotypeOneA.Concat(haplotypeTwoA).ToList();

        XYLinLen = YChrs.Select(c => (long)chrLengths[c]).Sum();
        XXLinLen = XChrs.Select(c => (long)chrLengths[c]).Sum();
        AutosomeLinLen = AutChrs.Select(c => (long)chrLengths[c]).Sum();
        XYGenomeLen = XYGenome.Sum(regions => regions.Sum(r => r.Length));
        XXGenomeLen = XXGenome.Sum(regions => regions.Sum(r => r.Length));
        AutosomeLen = Autosome.Sum(regions => regions.Sum(r => r.Length));

        GenContentsDict = genContentsDict;
    }
    
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

    public List<List<Region>> GetGenotype(SexType sexType)
        => sexType switch
        {
            SexType.Female => XXGenome,
            SexType.Male => XYGenome,
            _ => Autosome
        };
    private Dictionary<GeneListType, List<Gene> > GetRegionGenes(string chrNo)
        => GeneLists.ToDictionary(
            geneTypeList => geneTypeList.Key,
            geneTypeList => geneTypeList.Value
                    .Where(chromGenes => chromGenes.Key == chrNo)
                    .SelectMany(genes => genes.Value).ToList()
            );
    
    private Region GetRegion(string chrNo, bool isFirstHaplotype) 
        => new(0, ChrLengths[chrNo], chrNo, isFirstHaplotype,  new List<SNV>(), GetRegionGenes(chrNo));
    
    private IEnumerable<List<Region>> CreateHaplotype(SexType sex, bool firstHap)
        => ChrIDsForHap(sex, firstHap).Select(chr => new List<Region> { GetRegion(chr, firstHap) });
}