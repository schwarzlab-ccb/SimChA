using System.Collections.Immutable;
using System.Globalization;
using System.Text.Json;
using SimChA.Computation;
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
    private const string CONSISTENT_CNS_FILENAME = "consistent_CNs.tsv";
    private const string KARYOTYPES_FILENAME = "karyotypes.tsv";
    private const string CLONES_FILENAME = "clones.tsv";
    private const string CN_EVENTS_FILENAME = "events.tsv";
    private const string VCF_FILENAME = "vcf.tsv";
    // private const string FASTA_FILENAME = "genome.fa";
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
    
    public void WriteConsistentCNs(GenRef genRef, List<Sample> samples)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), CONSISTENT_CNS_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        
        outputFile.WriteLine(CopyNumbers.Header(true));
        var karyotypes = samples.Select(s => s.Karyotype);
        var breaks = karyotypes.Select(k => CopyNumbers.GetSegPoints(genRef, k)).ToList();
        var segmentation = genRef.AllChrs.ToDictionary(
            chrom => chrom, 
            chrom => breaks.SelectMany(br => br[chrom]).ToHashSet().OrderBy(val => val).ToList());
        foreach (var sample in samples)
        {
            var cns = CopyNumbers.CalcCNs(genRef, sample.Karyotype, segmentation);
            outputFile.WriteLine(CopyNumbers.ToTSV(cns, sample.SampleId));
        }
    }
    
    public void WriteCopyNumbers(GenRef genRef, List<Sample> samples)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), CN_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        outputFile.WriteLine(CopyNumbers.Header(true));

        foreach (var sample in samples)
        {
            var cns = CopyNumbers.CalcCNs(genRef, sample.Karyotype);
            outputFile.WriteLine(CopyNumbers.ToTSV(cns, sample.SampleId));
        }
    }

    public void WriteSimParams(SimChAConfig simChAConfig, string? name = null)
    {
        string filePath = (name != null) ? Path.Combine(Path.GetFullPath(OutFolder), name)
                                         : Path.Combine(Path.GetFullPath(OutFolder), SIM_PARAMS_FILENAME);
        using var file = new StreamWriter(filePath);
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
        // TODO: Do we need the header information?
        outputFile.WriteLine("##fileformat=VCFv4.3");
        outputFile.WriteLine("##source=SimChAV1.0");
        outputFile.WriteLine($"##reference=verily_{genRef.Name}_genome.fa");
        outputFile.WriteLine("#SAMPLEID\tCHROM\tPOS\tID\tREF\tALT");
        foreach (var sample in samples)
        {
            foreach (var snv in sample.Karyotype.GetFinalSNVs())
            {
                if (genRef.GenContentsDict == null)
                {
                    throw new Exception("Genomic Content hasn't been set correctly to allow SNV list to be created");
                }
                char refBase = genRef.GenContentsDict[snv.ChrNo][(int)snv.Location];
                // The VCF should *not* be aware of SNVs that didn't end up altering the location in the final karyotype
                if (char.ToUpper(refBase) != snv.Alt.ToString()[0])
                {
                    outputFile.WriteLine($"{sample.SampleId}\t{snv.ChrNo}\t{snv.Location}\t.\t{refBase}\t{snv.Alt}");
                }
            }
        }
    }

    public void WriteFasta(GenRef genRef, List<Sample> samples)
    {
        if (genRef.GenContentsDict == null)
        {
            throw new Exception("Reference Genome was not set. Please check that you have downloaded the correct assembly (see DownloadRefData.sh)");
        }
        foreach (var sample in samples)
        {
            string outPath = Path.Combine(Path.GetFullPath(OutFolder), $"{sample.SampleId}_genome.fa");
            Console.WriteLine($"Writing to file {outPath}");
            using var outputFile = new StreamWriter(outPath);
            var kar = sample.Karyotype;

            foreach (int contigId in kar.ContigIds())
            {
                outputFile.WriteLine($">ctg{contigId}");
                Console.WriteLine($"Writing out contig {contigId}");
                foreach (var region in kar.GetContig(contigId).GetRegions())
                {
                    string chrNo = region.ChrNo;
                    long start = region.Start;
                    long end   = region.End;
                    var regionSeq = new StringBuilder (genRef.GenContentsDict[chrNo].ToString((int)start, (int)(end-start)));
                    if (region.SNVs != null)
                    {
                        foreach (var snv in region.SNVs)
                        {
                            long loc = snv.Location - start;
                            regionSeq[(int)loc] = snv.Alt.ToString()[0];
                        }
                    }
                    if (!region.Forward)
                    {
                        char[] baseArray = regionSeq.ToString().ToCharArray();
                        Array.Reverse(baseArray);
                        regionSeq = new StringBuilder(new string(baseArray));
                    }
                    outputFile.Write(regionSeq);
                    // Note that the output FASTA will not be in blocks of 80 characters wide
                    // Implementing it means the write-out function takes way longer.
                }
                outputFile.Write("\n");
            }
        }
    }

    public void WriteSamples(List<SampleStat> cloneStats)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), CLONES_FILENAME);
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
        string fileFormat = filePath.Substring(filePath.Length - 3);
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

    // TODO: Is needed?
    public static List<(double fitness, int eventCount)> ReadFitnesses(string filePath, FitParams fitParams)
    {
        string fileFullPath = Path.GetFullPath(filePath);
        if (!File.Exists(fileFullPath))
        {
            throw new Exception($"File {fileFullPath} does not exist");
        }
        try
        {
            var fitnessFile = new StreamReader(fileFullPath);
            return Parsers.ParseClones(fitnessFile, fitParams);
        }
        catch (Exception e)
        {
            throw new Exception($"Failed to parse the file {fileFullPath}. Error {e.Message}");
        }
    }

    // TODO: Is needed?
    public static Dictionary<string, (double, double, double, int)> ReadCloneComponents(string filePath)
    {
        string fileFullPath = Path.GetFullPath(filePath);
        if (!File.Exists(fileFullPath))
        {
            throw new Exception($"File {fileFullPath} does not exist");
        }
        try
        {
            var fitnessFile = new StreamReader(fileFullPath);
            return Parsers.ParseCloneComponents(fitnessFile);
        }
        catch (Exception e)
        {
            throw new Exception($"Failed to parse the file {fileFullPath}. Error {e.Message}");
        }
    }

    public static Dictionary<string, StringBuilder> ReadFasta(List<string> allChrs, string folder)
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
    
    public static Dictionary<GeneListType, Dictionary<string, List<Gene>>> ReadGeneLists(string folder, Dictionary<string, SexType> chrSex)
    {
        var geneLists = new Dictionary<GeneListType, Dictionary<string, List<Gene>>>();
        var fileMap = new Dictionary<GeneListType, string>
        {
            { GeneListType.TumorSuppressor, TSGS_TSV },
            { GeneListType.Oncogene, OGS_TSV },
            { GeneListType.Essentiality, ESSENTIALS_TSV }
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
                geneLists[key] = Parsers.ParseGeneList(geneFile, chrNames);
            }
            catch (Exception e)
            {
                throw new Exception($"Failed to parse the file {fileFullPath}. Error {e.Message}");
            }
        }
        return geneLists;
    }

    public static  List<Sample> ReadProfiles(GenRef genRef, string cnaProfile, bool autosomesOnly)
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

    public static Dictionary<string, int> ReadEventCounts(string filePath)
    {
        string fileFullPath = Path.GetFullPath(filePath);
        if (!File.Exists(fileFullPath))
        {
            throw new Exception($"File {fileFullPath} does not exist");
        }
        try
        {
            var fitnessFile = new StreamReader(fileFullPath);
            return Parsers.ParseEventCounts(fitnessFile);
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

    public static GenRef GetGenRef(string dataFolder, bool useVariants = false)
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
