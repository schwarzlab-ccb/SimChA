using System.Globalization;
using System.Text.Json;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.Simulation;
using System.Text;

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
    private const string COPYNUMBERS_FILENAME = "copynumbers.tsv";
    private const string BINNED_COPYNUMBERS_FILENAME = "binned_CNs.tsv";
    private const string CONSISTENT_CNS_FILENAME = "consistent_CNs.tsv";
    private const string KARYOTYPES_FILENAME = "karyotypes.tsv";
    private const string CLONES_FILENAME = "clones.tsv";
    private const string CN_EVENTS_FILENAME = "events.tsv";
    private const string VCF_FILENAME = "vcf.tsv";
    private const string FASTA_FILENAME = "genome.fa";
    private const string SCORES_FILENAME = "optimization_scores.tsv";

    
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

    public void WriteSamples(IEnumerable<Sample> samples)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), SAMPLES_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        outputFile.WriteLine(Sample.Header());
        foreach (var sample in samples)
        {
            outputFile.WriteLine(sample.ToTSV());
        }
    }
    
    public void WriteConsistentCNs(GenRef genRef, IList<Sample> samples)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), CONSISTENT_CNS_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        
        outputFile.WriteLine("sample_id\tchrom\tstart\tend\tcn_a\tcn_b\tn_snvs");

        var segs = CopyNumbers.GetSegPoints(genRef , samples.SelectMany(s => s.Kars.Values).ToList());
        
        foreach (var sample in samples)
        {
            foreach (var clone in sample.Clones)
            {
                var kar = sample.Kars[clone.CloneId];
                var cns = CopyNumbers.CalcCopyNumbers(genRef, kar, segs, sample.SexXX, true);
                string name = sample.Clones.Count > 1 ? $"{sample.SampleId}_{clone.CloneId}" : $"{sample.SampleId}";
                outputFile.WriteLine(CopyNumbers.ToTSV(cns, name, false));
            }
        }
    }
    
    public void WriteCopyNumbers(GenRef genRef, IEnumerable<Sample> samples)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), COPYNUMBERS_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        outputFile.WriteLine("sample_id\tchrom\tstart\tend\tcn_a\tcn_b\tn_snvs");

        foreach (var sample in samples)
        {
            foreach (var clone in sample.Clones)
            {
                var cns = CopyNumbers.CalcCopyNumbers(genRef, sample.Kars[clone.CloneId], sample.SexXX);
                string name = sample.Clones.Count > 1 ? $"{sample.SampleId}_{clone.CloneId}" : $"{sample.SampleId}";
                outputFile.WriteLine(CopyNumbers.ToTSV(cns, name, false));
            }
        }
    }

    public void WriteCopyNumbers(Dictionary<string, List<CopyNumber>> cnProfiles)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), BINNED_COPYNUMBERS_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        outputFile.WriteLine("sample_id\tchrom\tstart\tend\tcn_a\tcn_b\tn_snvs");
        foreach (var cnProfile in cnProfiles)
        {
            outputFile.WriteLine(CopyNumbers.ToTSV(cnProfile.Value, cnProfile.Key, false));
        }
    }

    public void WriteSimParams(SimParams simParams, string? name = null)
    {
        string filePath = (name != null) ? Path.Combine(Path.GetFullPath(OutFolder), name)
                                         : Path.Combine(Path.GetFullPath(OutFolder), SIM_PARAMS_FILENAME);
        using var file = new StreamWriter(filePath);
        var options = new JsonSerializerOptions { IncludeFields = true, WriteIndented = true };
        string jsonString = JsonSerializer.Serialize(simParams, options);
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
            foreach (var kar in sample.Kars)
            {
                string sampleName = sample.Clones.Count > 1 ? $"{sample.SampleId}_{kar.Key}" : $"{sample.SampleId}";
                outputFile.WriteLine($"{sampleName}\t{kar.Value}");
            }
        }
    }

    public void WriteEvents(IEnumerable<Sample> samples)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), CN_EVENTS_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        outputFile.WriteLine("sample_id\tevent_type\tdepth\tevent_string\tdelta_fitness\ttotal_fitness");
        foreach (var sample in samples)
        {
            foreach (var clone in sample.EventDescs)
            {
                foreach (var cnEvent in clone.Value)
                {   
                    string sampleName = sample.Clones.Count > 1 ? $"{sample.SampleId}_{clone.Key}" : $"{sample.SampleId}";
                    outputFile.WriteLine($"{sampleName}\t{cnEvent.EventType}\t{cnEvent.Depth}\t{cnEvent.Description}" +
                                         $"\t{cnEvent.DeltaFitness:f6}\t{cnEvent.TotalFitness:f6}");
                }
            }
        }
    }

    public void WriteScores(IEnumerable<Dictionary<string, double>> scores)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), SCORES_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        var keyOrder = scores.First().Keys.ToList();
        var header = string.Join("\t", keyOrder);
        outputFile.WriteLine(header);
        foreach (var score in scores)
        {
            string scoreline = "";
            foreach (var key in keyOrder)
            {
                scoreline += $"{score[key]}\t";
            }
            scoreline = scoreline.Substring(0, scoreline.Length - 2);
            outputFile.WriteLine(scoreline);
        }
    }

    public void WriteFitnessLandscape(string filename, List<List<double>> output)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), filename);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        outputFile.WriteLine("alpha\tbeta\tfitness");
        foreach (var line in output)
        {
            outputFile.WriteLine($"{line[0]}\t{line[1]}\t{line[2]}");
        }
    }

    public void WriteVCF(GenRef genRef, IEnumerable<Sample> samples)
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
            foreach (var clone in sample.EventDescs)
            {
                var kar = sample.Kars[clone.Key];
                string sampleName = sample.Clones.Count > 1 ? $"{sample.SampleId}_{clone.Key}" : $"{sample.SampleId}";
                var finalSNVs = kar.GetFinalSNVs();

                foreach (var snv in finalSNVs)
                {
                    if (genRef.GenContentsDict == null)
                    {
                        throw new Exception("Genomic Content hasn't been set correctly to allow SNV list to be created");
                    }
                    char refBase = genRef.GenContentsDict[snv.chrNo][(int)snv.location];
                    // The VCF should *not* be aware of SNVs that didn't end up altering the location in the final karyotype
                    if (char.ToUpper(refBase) != snv.newBase.ToString()[0])
                    {
                        outputFile.WriteLine($"{sampleName}\t{snv.chrNo}\t{snv.location}\t.\t{refBase}\t{snv.newBase}");
                    }
                }
            }
        }
    }

    public void WriteFasta(GenRef genRef, IEnumerable<Sample> samples)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), FASTA_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        // TODO: Do we want WriteFasta to work with multiple samples? Currently only set up for single samples
        var count = 0;
        if (genRef.GenContentsDict == null)
        {
            throw new Exception("Reference Genome was not set. Please check that you have downloaded the correct assembly (see DownloadRefData.sh)");
        }
        foreach (var sample in samples)
        {
            foreach (var clone in sample.EventDescs)
            {
                var kar = sample.Kars[clone.Key];

                foreach (var contigId in kar.ContigIds())
                {
                    outputFile.WriteLine($">ctg{contigId}");
                    Console.WriteLine($"Writing out contig {contigId}");
                    foreach (var region in kar.GetContig(contigId).GetRegions())
                    {
                        var chrNo = region.ChrNo;
                        var start = region.Start;
                        var end   = region.End;
                        var regionSeq = new StringBuilder (genRef.GenContentsDict[chrNo].ToString((int)start, (int)(end-start)));
                        if (region.SNVDict != null)
                        {
                            foreach (var snv in region.SNVDict)
                            {
                                var loc = snv.Key - start;
                                regionSeq[(int)(snv.Key-start)] = snv.Value.ToString()[0];
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
                if (count > 0) return;
            }
        }
    }

    public void WriteClones(IEnumerable<Sample> samples)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), CLONES_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var file = new StreamWriter(outPath);
        file.WriteLine("sample_id\tploidy\tcoverage\tfitness\tstress\ttsg\tog\tess");
        foreach (var sample in samples)
        {
            foreach (var stats in sample.Stats)
            {
                string sampleName = sample.Clones.Count > 1 ? $"{sample.SampleId}_{stats.Key}" : $"{sample.SampleId}";
                var clone = stats.Value;
                file.WriteLine($"{sampleName}\t{clone.Ploidy}\t{clone.Coverage}\t{clone.Fitness}\t{clone.Stress}\t{clone.Tsg}\t{clone.Og}\t{clone.Ess}");
            }
        }
    }

    public static List<CloneIn> ReadClones(string filePath, bool parseFitness)
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
            return Parsers.ParseClones(cloneFile, parseFitness, separator);
        }
        catch (Exception e)
        {
            throw new Exception($"Failed to parse the file {fileFullPath}. Error {e.Message}");
        }
    }

    public static List<(double fitness, int eventCount)> ReadFitnesses(string filePath, FitnessParams fitnessParams)
    {
        string fileFullPath = Path.GetFullPath(filePath);
        if (!File.Exists(fileFullPath))
        {
            throw new Exception($"File {fileFullPath} does not exist");
        }
        try
        {
            var fitnessFile = new StreamReader(fileFullPath);
            return Parsers.ParseClones(fitnessFile, fitnessParams);
        }
        catch (Exception e)
        {
            throw new Exception($"Failed to parse the file {fileFullPath}. Error {e.Message}");
        }
    }

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

    public static SimParams ReadSimParams(string filePath)
    {
        string fileFullPath = Path.GetFullPath(filePath);
        if (!File.Exists(fileFullPath))
        {
            throw new Exception($"File {fileFullPath} does not exist");
        }
        try
        {
            string serializedJSON = File.ReadAllText(fileFullPath);
            return Parsers.ParseSimParams(serializedJSON);
        }
        catch (Exception e)
        {
            throw new Exception($"Failed to parse the file {fileFullPath}. Error {e.Message}");
        }
    }
    
    public static Dictionary<GeneListType, Dictionary<string, List<Gene>>> ReadGeneLists(string folder, Dictionary<string, SexEnum> chrSex)
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

    public static Dictionary<string, Karyotype> ReadProfiles(GenRef genRef, string cnaProfile)
    {
        string fileFullPath = Path.GetFullPath(cnaProfile);
        if (!File.Exists(fileFullPath))
        {
            throw new Exception($"File {fileFullPath} does not exist");
        }
        try
        {
            var cnaFile = new StreamReader(fileFullPath);
            var profiles = Parsers.ParseCNAProfile(genRef, cnaFile);
            foreach (var pro in profiles)
            {
                pro.Value.GlueNeighbours();
            }
            return profiles;
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

    public static Dictionary<string, List<CopyNumber>> ReadProfiles(string cnaProfile)
    {
        string fileFullPath = Path.GetFullPath(cnaProfile);
        if (!File.Exists(fileFullPath))
        {
            throw new Exception($"File {fileFullPath} does not exist");
        }
        try
        {
            var cnaFile = new StreamReader(fileFullPath);
            return Parsers.ParseCNAProfile(cnaFile);
        }
        catch (Exception e)
        {
            throw new Exception($"Failed to parse the file {fileFullPath}. Error {e.Message}");
        }
    }
    
    public static (Dictionary<string, int> chrLengths, Dictionary<string, SexEnum> chrSex) ReadChromosomes(string folder)
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

    public static (Dictionary<string, (int, int)> p, Dictionary<string, (int, int)> q) ReadCentromeres(string folder)
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

    public static GenRef GetGenRef(string dataFolder, bool includeSexChromosomes = true, bool useVariants = false)
    {
        string refName = Path.GetFileName(dataFolder);
        var (chrLengths, chrSex)  = ReadChromosomes(dataFolder);
        var centromeres = ReadCentromeres(dataFolder);
        var allChrs = chrSex.Select(pair => pair.Key).ToList();
        var genContentsDict = useVariants ? ReadFasta(allChrs, dataFolder) : null;
        var geneLists = ReadGeneLists(dataFolder, chrSex);
        return new GenRef(refName, chrLengths, chrSex, centromeres, geneLists, includeSexChromosomes, genContentsDict);
    }
}
