using System.Globalization;
using System.Text;
using System.Text.Json;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;
using SimChA.Simulation;

namespace SimChA.IO;

public class FileIO
{
    private const string NEWICK_FILENAME = "parent_graph.new";
    private const string SAMPLES_FILENAME = "samples.tsv";
    private const string SUBCLONES_FILENAME = "clones.tsv";
    private const string COPYNUMBERS_FILENAME = "copynumbers.tsv";
    private const string SIM_PARAMS_FILENAME = "sim_params.json";
    private const string CN_EVENTS_FILENAME = "events.tsv";
    private const string ESSENTIALS_TSV = "essentials.tsv";
    private const string OGS_TSV = "OGs.tsv";
    private const string TSGS_TSV = "TSGs.tsv";
    private const string SAMPLE_FITNESS_FILE = "sample_fitness.tsv";
    private const string FITNESS_INPUT_FILE = "subclones.out";
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
            outputFile.WriteLine(sample.ToString());
        }
    }

    public void WriteClones(IEnumerable<Sample> samples)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), SUBCLONES_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);

        outputFile.WriteLine("sample_id\t" + CloneIn.Header());
        foreach (var sample in samples)
        {
            foreach (var clone in sample.Clones)
            {
                outputFile.WriteLine($"{sample.SampleId}\t{clone}");
            }
        }
    }

    public void WriteCopyNumbers(IEnumerable<Sample> samples)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), COPYNUMBERS_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);

        var copyNumbersString = new StringBuilder();
        copyNumbersString.Append("sample_id\tclone_id\tchr\tstart\tend\tcn_a\tcn_b\n");

        foreach (var sample in samples)
        {
            foreach (var clone in sample.Clones)
            {
                var copynumbers = CopyNumbers.CalcCopyNumbers(sample.Kars[clone.CloneId], sample.SexXX);
                string name = sample.Clones.Count > 1 ? $"{sample.SampleId}_{clone.CloneId}" : $"{sample.SampleId}";
                copyNumbersString.Append(CopyNumbers.ToTSV(copynumbers, name, false) + "\n");
            }
            outputFile.Write(copyNumbersString.ToString());
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

    public void WriteEvents(IEnumerable<Sample> samples)
    {
        //TODO: Format output, talk with Tom about readable ideas
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), CN_EVENTS_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        StringBuilder abberationString = new();
        abberationString.Append("sample\tid\tevent_type\tevent_string\tdelta_fitness\ttotal_Fitness\tmutations\n");
        foreach (var sample in samples)
        {
            foreach (var cnEvent in sample.EventDescs)
            {
                abberationString.Append($"{cnEvent.CloneId}\t" +
                                        $"{cnEvent.EventType}\t" +
                                        $"{cnEvent.Description}\t" +
                                        $"{Math.Round((decimal)cnEvent.DeltaFitness, 8)}\t" +
                                        $"{Math.Round((decimal)cnEvent.TotalFitness, 8)}\t");
            }
        }
        outputFile.Write(abberationString.ToString());
    }

    public void WriteFitness(Dictionary<string, List<CloneStat>> sampleStats)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), SAMPLE_FITNESS_FILE);
        Console.WriteLine($"Writing to file {outPath}");
        using var file = new StreamWriter(outPath);
        var myT = typeof(CloneStat);
        var fileds = myT.GetProperties();
        var fieldNames = fileds.Select(f => f.Name).ToList();
        file.WriteLine("sample_id\t" + string.Join("\t", fieldNames));
        foreach (var sample in sampleStats)
        {
            foreach (var clone in sample.Value)
            {
                // Get all the field values of the record sample
                var values = fileds.Select(f => $"{f.GetValue(clone):f4}");
                file.WriteLine(sample.Key + "\t" + string.Join("\t", values));
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

    public static string ReadNewick(string newickFile)
    {
        string fileFullPath = Path.GetFullPath(newickFile);
        if (!File.Exists(fileFullPath))
        {
            throw new Exception($"File {fileFullPath} does not exist");
        }
        try
        {
            return File.ReadAllText(fileFullPath);
        }
        catch (Exception e)
        {
            throw new Exception($"Failed to parse the file {fileFullPath}. Error {e.Message}");
        }
    }

    public static Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> ReadGeneLists(
        string folder,
        bool isFemale,
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
                geneLists[key] = Parsers.ParseGeneList(geneFile, isFemale);
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