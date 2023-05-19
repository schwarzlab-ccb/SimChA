using System.Globalization;
using System.Text;
using System.Text.Json;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.IO;

public class FileIO
{
    private const string SAMPLES_FILENAME = "samples.tsv";
    private const string COPYNUMBERS_FILENAME = "copynumbers.tsv";
    private const string SIM_PARAMS_FILENAME = "sim_params.json";
    private const string KARYOTYPES_FILENAME = "karyotypes.tsv";
    private const string CLONES_FILENAME = "clones.tsv";
    private const string CN_EVENTS_FILENAME = "events.tsv";
    private const string ESSENTIALS_TSV = "essentials.tsv";
    private const string OGS_TSV = "OGs.tsv";
    private const string TSGS_TSV = "TSGs.tsv";
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
        else
        {
            foreach (string file in Directory.EnumerateFiles(OutFolder))
            {
                File.Delete(file);
            }
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
    
    public void WriteCopyNumbers(IEnumerable<Sample> samples)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), COPYNUMBERS_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        
        outputFile.WriteLine("sample_id\tchrom\tstart\tend\tcn_a\tcn_b");

        foreach (var sample in samples)
        {
            foreach (var clone in sample.Clones)
            {
                var cns = CopyNumbers.CalcCopyNumbers(sample.Kars[clone.CloneId], sample.SexXX);
                string name = sample.Clones.Count > 1 ? $"{sample.SampleId}_{clone.CloneId}" : $"{sample.SampleId}";
                outputFile.WriteLine(CopyNumbers.ToTSV(cns, name, false));
            }
        }
    }

    public void WriteSimParams(SimParams simParams)
    {
        string filePath = Path.Combine(Path.GetFullPath(OutFolder), SIM_PARAMS_FILENAME);
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
        // TODO: Format output, talk with Tom about readable ideas
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
        if (!File.Exists(fileFullPath))
        {
            throw new Exception($"File {fileFullPath} does not exist");
        }
        try
        {
            var cloneFile = new StreamReader(fileFullPath);
            return Parsers.ParseClones(cloneFile, parseFitness);
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
    
    public static Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> ReadGeneLists(
        string folder,
        GenomeAssembly assembly)
    {
        var geneLists = new Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>>();
        var fileMap = new Dictionary<GeneListType, string>
        {
            { GeneListType.TumorSuppressor, TSGS_TSV },
            { GeneListType.Oncogene, OGS_TSV },
            { GeneListType.Essentiality, ESSENTIALS_TSV }
        };
        foreach ((var key, string filename) in fileMap)
        {
            string filePath = Path.Combine(folder, assembly.ToString(), filename);
            string fileFullPath = Path.GetFullPath(filePath);
            if (!File.Exists(fileFullPath))
            {
                throw new Exception($"File {fileFullPath} does not exist");
            }
            try
            {
                var geneFile = new StreamReader(fileFullPath);
                geneLists[key] = Parsers.ParseGeneList(geneFile);
            }
            catch (Exception e)
            {
                throw new Exception($"Failed to parse the file {fileFullPath}. Error {e.Message}");
            }
        }
        return geneLists;
    }

    public static Dictionary<string, Karyotype> ReadProfiles(string cnaProfile)
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
}