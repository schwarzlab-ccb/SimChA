using System.Globalization;
using System.Text.Json;
using System.Text;
using SimChA.Data;
using SimChA.EventData;

namespace SimChA.IO;

public class FileIO
{
    // data
    private const string CHROMOSOMES_TSV = "chromosomes.tsv";
    private const string ESSENTIALS_TSV = "essentials_select.tsv";
    private const string OGS_TSV = "ogs_select.tsv";
    private const string TSGS_TSV = "tsgs_select.tsv";
    private const string GENOME_FASTA = "genome.fa";
    private const string CENTROMERES_TSV = "centromeres.tsv";
    // input
    private const string SIM_PARAMS_FILENAME = "sim_params.json";
    
    // output
    private const string SAMPLES_FILENAME = "samples.tsv";
    private const string CN_FILENAME = "copynumbers.tsv";
    private const string KARYOTYPES_FILENAME = "karyotypes.tsv";
    private const string CN_EVENTS_FILENAME = "events.tsv";
    private const string VCF_FILENAME = "vcf.tsv";
    // private const string FITNESSES_FILENAME = "mcmc_fitnesses.tsv";
    // private const string TREE_FILENAME = "tree.tsv";
    
    private string Timestamp { get; }
    private string OutFolder { get; }

    public FileIO(string outFolder)
    {
        Timestamp = DateTime.Now.ToString("yy_MM_dd_HH_mm_ss");
        CultureInfo.DefaultThreadCurrentCulture = CultureInfo.InvariantCulture;

        OutFolder = outFolder;
        if (!Directory.Exists(OutFolder))
        {
            Directory.CreateDirectory(OutFolder);
        }
    }
    
    public void WriteCopyNumbers(Dictionary<string, IEnumerable<CopyNumber>> samples)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), CN_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        outputFile.WriteLine("sample_id\t" + CopyNumber.Header());
        foreach ((string sampleID, var cns) in samples)
        {
            foreach (var cn in cns)
            {
                outputFile.WriteLine($"{sampleID}\t{cn.ToTSV()}");
            }
        }
    }

    public void WriteSimParams(SimChAConfig simChAConfig)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), SIM_PARAMS_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var file = new StreamWriter(outPath);
        var options = new JsonSerializerOptions { IncludeFields = true, WriteIndented = true };
        string jsonString = JsonSerializer.Serialize(simChAConfig, options);
        file.WriteLine(jsonString);
    }

    public void WriteKaryotypes(List<Sample> samples)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), KARYOTYPES_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        outputFile.WriteLine("sample_id\tkaryotype");
        
        foreach (var sample in samples)
        {
            var kar = sample.Karyotype;
            outputFile.WriteLine($"{sample.SampleId}\t{kar}");
        }
    }

    public void WriteEvents(List<Sample> samples)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), CN_EVENTS_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        outputFile.WriteLine("sample_id\t" + CNEventDesc.Header());
        foreach (var sample in samples)
        {
            foreach (var cnEvent in sample.Events)
            {   
                outputFile.WriteLine($"{sample.SampleId}\t{cnEvent.ToTSV()}");
            }
        }
    }
    
    public void WriteVCF(GenRef genRef, List<Sample> samples)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), VCF_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        outputFile.WriteLine("##fileformat=VCFv4.3");
        outputFile.WriteLine("##source=SimChAV1.0");
        outputFile.WriteLine($"##reference=verily_{genRef.Name}_genome.fa");
        outputFile.WriteLine("#SAMPLEID\tCHROM\tPOS\tID\tREF\tALT");
        foreach (var sample in samples)
        {
            foreach (var snv in sample.Karyotype.GetSNVs())
            {
                var refBase = genRef.GetRefBase(snv.Chrom, (int)snv.Pos);
                // The VCF should *not* be aware of SNVs that didn't end up altering the location in the final karyotype
                if (refBase != snv.Alt)
                {
                    outputFile.WriteLine($"{sample.SampleId}\t{snv.Chrom}\t{snv.Pos}\t.\t{refBase}\t{snv.Alt}");
                }
            }
        }
    }

    public void WriteFasta(List<Sample> samples)
    {
        foreach (var sample in samples)
        {
            string outPath = Path.Combine(Path.GetFullPath(OutFolder), $"{sample.SampleId}_genome.fa");
            Console.WriteLine($"Writing to file {outPath}");
            using var outputFile = new StreamWriter(outPath);
            foreach (string region in sample.Karyotype.GetSeq())
            {
                outputFile.Write(region);
            }
        }
    }

    public void WriteSamples(List<SampleStat> cloneStats)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), SAMPLES_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var file = new StreamWriter(outPath);
        file.WriteLine(SampleStat.Header());

        foreach (var clone in cloneStats)
        {
            file.WriteLine(clone.ToString());
        }
    }

    public static (CTreeNode root, List<CTreeNode> tree) ReadCloneTree(string filePath, bool parseFitness)
    {
        string fileFullPath = Path.GetFullPath(filePath);
        string fileFormat = filePath[^3..];
        if (fileFormat != "tsv" && fileFormat != "csv")
        {
            throw new Exception($"File {filePath} should be a tsv or csv.");
        }
        string separator = fileFormat == "tsv" ? "\t" : ",";
        if (!File.Exists(fileFullPath))
        {
            throw new Exception($"File {fileFullPath} does not exist");
        }
        try
        {
            var cloneFile = new StreamReader(fileFullPath);
            var tree = Parsers.ParseClonesWithEvents(cloneFile, parseFitness, separator);
            var selfParent = tree.FindAll(n => n.ParentId == n.CloneId);
            return selfParent.Count switch
            {
                > 1 => throw new Exception($"More than one ({selfParent.Count}) root nodes (parented to self) found in the clone tree {filePath}."),
                0 => throw new Exception($"No root node found in the clone tree  {filePath}."),
                _ => (selfParent[0], tree)
            };
        }
        catch (Exception e)
        {
            throw new Exception($"Failed to parse the file {fileFullPath}. Error {e.Message}");
        }
    }

    private static Dictionary<string, StringBuilder> ReadFasta(List<string> allChrs, string folder)
    {
        string fileFullPath = Path.GetFullPath(Path.Combine(folder, GENOME_FASTA));
        if (!File.Exists(fileFullPath))
        {
            throw new Exception($"File {fileFullPath} does not exist");
        }
        try
        {
            var fastaFile = new StreamReader(fileFullPath);
            var sequenceList =  Parsers.ParseFasta(fastaFile).ToList();
            return sequenceList.Select((seq, i) => new {i, seq}).ToDictionary(x => allChrs[x.i], x => x.seq);
        }
        catch (Exception e)
        {
            throw new Exception($"Failed to parse the file {fileFullPath}. Error {e.Message}");
        }
    }

    public static SimChAConfig ReadSimChAConfig(string filePath)
    {
        string fileFullPath = Path.GetFullPath(filePath);
        if (!File.Exists(fileFullPath))
        {
            throw new Exception($"File {fileFullPath} does not exist");
        }
        try
        {
            string serializedJSON = File.ReadAllText(fileFullPath);
            var config = Parsers.ParseSimParams(serializedJSON);
            if (config.SimParams == null)
            {
                throw new Exception("No simulation parameters found. Please set \"SimParams\" in the config JSON.");
            }
            if (config.FitParams == null)
            {
                throw new Exception("No fitness parameters found. Please set \"FitParams\" in the config JSON.");
            }
            return config;
        }
        catch (Exception e)
        {
            throw new Exception($"Failed to parse the file {fileFullPath}. Error {e.Message}");
        }
    }

    private static Dictionary<GeneLT, Dictionary<string, List<Gene>>> ReadGeneLists(string folder, Dictionary<string, SexType> chrSex)
    {
        var geneLists = new Dictionary<GeneLT, Dictionary<string, List<Gene>>>();
        var fileMap = new Dictionary<GeneLT, string>
        {
            { GeneLT.TSG, TSGS_TSV },
            { GeneLT.OG, OGS_TSV },
            { GeneLT.Ess, ESSENTIALS_TSV }
        };
        foreach ((var key, string filename) in fileMap)
        {
            string filePath = Path.Combine(folder, filename);
            string fileFullPath = Path.GetFullPath(filePath);
            if (!File.Exists(fileFullPath))
            {
                throw new Exception($"File {fileFullPath} does not exist");
            }
            try
            {
                var geneFile = new StreamReader(fileFullPath);
                var chrNames = chrSex.Select(pair => pair.Key).ToList();
                geneLists[key] = Parsers.ParseGeneList(geneFile, chrNames, key);
            }
            catch (Exception e)
            {
                throw new Exception($"Failed to parse the file {fileFullPath}. Error {e.Message}");
            }
        }
        return geneLists;
    }

    public static List<Sample> ReadProfiles(GenRef genRef, string cnaProfile, bool autosomesOnly)
    {
        string fileFullPath = Path.GetFullPath(cnaProfile);
        if (!File.Exists(fileFullPath))
        {
            throw new Exception($"File {fileFullPath} does not exist");
        }
        try
        {
            var cnaFile = new StreamReader(fileFullPath);
            var profiles = Parsers.ParseCNAProfile(genRef, cnaFile, autosomesOnly);
            var samples = new List<Sample>();
            foreach ((string sampleId, var karyotype) in profiles)
            {
                samples.Add(new Sample(sampleId, sampleId, karyotype));
            }
            return samples;
        }
        catch (Exception e)
        {
            throw new Exception($"Failed to parse the file {fileFullPath}. Error {e.Message}");
        }
    }
    
    private static (Dictionary<string, int> chrLengths, Dictionary<string, SexType> chrSex) ReadChromosomes(string folder)
    {
        string fileFullPath = Path.GetFullPath(Path.Combine(folder, CHROMOSOMES_TSV));
        if (!File.Exists(fileFullPath))
        {
            throw new Exception($"File {fileFullPath} does not exist");
        }
        try
        {
            string fileContent = File.ReadAllText(fileFullPath);
            return Parsers.ParseChromosomes(fileContent);
        }
        catch (Exception e)
        {
            throw new Exception($"Failed to parse the file {fileFullPath}. Error {e.Message}");
        }
    }

    private static Dictionary<string, GenRange> ReadCentromeres(string folder)
    {
        string fileFullPath = Path.GetFullPath(Path.Combine(folder, CENTROMERES_TSV));
        if (!File.Exists(fileFullPath))
        {
            throw new Exception($"File {fileFullPath} does not exist");
        }
        try
        {
            var fileContent = new StreamReader(fileFullPath);
            return Parsers.ParseCentromeres(fileContent);
        }
        catch (Exception e)
        {
            throw new Exception($"Failed to parse the file {fileFullPath}. Error {e.Message}");
        }
    }

    public static GenRef ReadGenRef(string dataFolder, bool useVariants = false)
    {
        string refName = Path.GetFileName(dataFolder);
        var (chrLengths, chrSex) = ReadChromosomes(dataFolder);
        var centromeres = ReadCentromeres(dataFolder);
        var allChrs = chrSex.Select(pair => pair.Key).ToList();
        var genContentsDict = useVariants ? ReadFasta(allChrs, dataFolder) : null;
        var geneLists = ReadGeneLists(dataFolder, chrSex);
        return new GenRef(refName, chrLengths, chrSex, centromeres, geneLists, genContentsDict);
    }
}
