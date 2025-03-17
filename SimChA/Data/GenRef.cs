using System.Text;

namespace SimChA.Data;

public class GenRef
{
    public string Name { get; }
    public Dictionary<string, int> ChrLengths { get; }
    private Dictionary<string, SexType> ChrSex { get; }
    public Dictionary<string, GenRange> Centromeres { get; }
    public string YChrName { get; }
    public string XChrName { get; }
    public List<List<List<Region>>> Genomes { get; }
    public List<long> GenomeLens { get; }
    public List<long> ReferenceLens { get; }
    public List<List<string>> SexChromNames { get; }
    public List<string> AllChrNames => SexChromNames[(int)SexType.Male];
    public int AutosomesCount => SexChromNames[(int) SexType.Any].Count;
    private List<Dictionary<string, List<Gene>>> GeneChromMap { get; }
    public List<List<List<Gene>>> GeneLists { get; }
    private Dictionary<string, StringBuilder>? GenContentsDict { get; }
    
    public GenRef(
        string name,
        Dictionary<string, int> chrLengths,
        Dictionary<string, SexType> chrSex,
        Dictionary<string, GenRange> centromeres,
        List<Dictionary<string, List<Gene>>> geneChromMap,
        Dictionary<string, StringBuilder>? genContentsDict = null)
    {
        Name = name;
        ChrLengths = chrLengths;
        ChrSex = chrSex;
        Centromeres = centromeres;
        GeneChromMap = geneChromMap;
        GenContentsDict = genContentsDict;

        YChrName = ChrSex.Where(pair => pair.Value == SexType.Male).Select(pair => pair.Key).FirstOrDefault("");
        XChrName = ChrSex.Where(pair => pair.Value == SexType.Female).Select(pair => pair.Key).FirstOrDefault("");
        SexChromNames = [
            ChrSex.Where(pair => pair.Value == SexType.Any).Select(pair => pair.Key).ToList(),
            ChrSex.Where(pair => pair.Value != SexType.Male).Select(pair => pair.Key).ToList(),
            ChrSex.Select(pair => pair.Key).ToList()
        ];
        Genomes = Enum.GetValues(typeof(SexType)).Cast<SexType>().Select(t
            => CreateHaplotype(t, true).Concat(CreateHaplotype(t, false)).ToList())
            .ToList();
        GenomeLens = Genomes.Select(genome 
            => genome.Sum(regs => regs.Sum(reg => reg.Length)))
            .ToList();
        ReferenceLens = SexChromNames.Select(sexList
            => sexList.Sum(chrName => (long) chrLengths[chrName]))
            .ToList();
        GeneLists = SexChromNames.Select(chrs => GeneChromMap.Select(map 
                => CreateGeneList(chrs, map)).ToList())
            .ToList();
    }

    private static List<Gene> CreateGeneList(List<string> chrs, Dictionary<string, List<Gene>> geneChromMap)
        => chrs.SelectMany(chrom => geneChromMap[chrom]).OrderBy(g => g.GeneId).ToList();

    public char[] GetGenContents(string chrom, int start, int length)
        => GenContentsDict == null
            ? ['N']
            : [..GenContentsDict[chrom].ToString(start, length)];

    private static Nucleotide CharToNucleotide(char c)
        => char.ToUpper(c) switch
        {
            'A' => Nucleotide.A,
            'C' => Nucleotide.C,
            'G' => Nucleotide.G,
            'T' => Nucleotide.T,
            _ => Nucleotide.N // Default case for invalid characters
        };

    public Nucleotide GetRefBase(string chrom, int location)
        => GenContentsDict == null
            ? Nucleotide.N
            : CharToNucleotide(GenContentsDict[chrom].ToString(location, 1)[0]);

    public Nucleotide GetRefBaseFromSeq(IEnumerable<string> seq, int location)
        => GenContentsDict == null 
            ? Nucleotide.N 
            : CharToNucleotide(seq.ElementAt(location)[0]);
    
    private List<string> ChrNamesForHap(SexType sex, bool firstHaplotype = true)
        => (firstHaplotype, sex) switch
        {
            (true, SexType.Male) => SexChromNames[(int) SexType.Female],
            (false, SexType.Male) => SexChromNames[(int) SexType.Any].Concat([YChrName]).ToList(),
            _ => SexChromNames[(int) sex]
        };

    private int GetExpCount(string chrom, SexType sexType)
    {
        if (chrom == XChrName)
        {
            return sexType switch
            {
                SexType.Female => 2,
                SexType.Male => 1,
                _ => 0
            };
        }
        if (chrom == YChrName)
        {
            return sexType == SexType.Male ? 1 : 0;
        }
        return 2;
    }
    
    public List<int[]> GetInitialGeneCounts(SexType sexType, bool keepEmpty)
    {
        List<int[]> res = [];
        foreach (var genesPerChrom in GeneChromMap)
        {
            int geneCount = genesPerChrom.Sum(pair => pair.Value.Count);
            int[] typeCounts = new int[geneCount];
            if (!keepEmpty)
            {
                foreach ((string chrom, var genes) in genesPerChrom)
                {
                    // In male the sex chromosomes have CN == 1, otherwise all are 2
                    int expCount = keepEmpty ? 0: GetExpCount(chrom, sexType);
                    foreach (var gene in genes)
                    {
                        typeCounts[gene.GeneId] = expCount;
                    }
                }
            }
            res.Add(typeCounts);
        }
        return res;
    }

    private IEnumerable<Gene> GetChromGenes(string chrNo)
        => GeneChromMap.SelectMany(chrGenes => chrGenes[chrNo]);

    private Region GetRegion(string chrNo, bool isFirstHaplotype)
        => new(0, ChrLengths[chrNo], chrNo, isFirstHaplotype, null, GetChromGenes(chrNo).ToList());

    private IEnumerable<List<Region>> CreateHaplotype(SexType sex, bool firstHap)
        => ChrNamesForHap(sex, firstHap).Select(chr => new List<Region> { GetRegion(chr, firstHap) });

    public IEnumerable<Gene> GetGenesBetween(string chrNo, int start, int end)
        => GeneChromMap.SelectMany(chrGenes => chrGenes[chrNo]
            .SkipWhile(g => g.Start < start)
            .TakeWhile(g => g.End <= end)
        );
}