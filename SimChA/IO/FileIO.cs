using System.Globalization;
using System.Text;
using System.Text.Json;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.IO;

public class FileIO
{
    private const string DOT_FILENAME = "parent_graph.dot";
    private const string NEWICK_FILENAME = "parent_graph.new";
    private const string SUBCLONES_FILENAME = "subclones.out";
    private const string COPYNUMBERS_FILENAME = "copynumbers.out";
    private const string BAF_FILENAME = "baf.out";
    private const string LOGR_FILENAME = "logr.out";
    private const string SIM_PARAMS_FILENAME = "sim_params.json";
    private const string CN_EVENTS_FILENAME = "abberations.tsv";
    private const string ESSENTIALS_TSV = "essentials.tsv";
    private const string OGS_TSV = "ogs.tsv";
    private const string TSGS_TSV = "tsgs.tsv";
    private const string SAMPLE_FITNESS_FILE = "sample_fitness.tsv";

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

    public void WriteClones(IEnumerable<Clone> subClones)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), SUBCLONES_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);

        var clonesString = new StringBuilder();
        foreach (var subClone in subClones)
        {
            clonesString.AppendLine(subClone.ToString());
        }

        outputFile.Write(clonesString.ToString());
    }
    
    public void WriteCopyNumbers(IEnumerable<Clone> subClones, bool isFemale)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), COPYNUMBERS_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);

        var copyNumbersString = new StringBuilder();
        copyNumbersString.Append("sample_id\tchr\tstart\tend\tcn_a\tcn_b\n");

        foreach (var subClone in subClones)
        {
            var copynumbers = CopyNumbers.CalcCopyNumbers(subClone.Karyotype, isFemale);
            copyNumbersString.Append(CopyNumbers.ToTSV(copynumbers, subClone.Name, false) + "\n");
        }
        outputFile.Write(copyNumbersString.ToString());
    }

    public void WriteSimParams(SimParams simParams)
    {
        string filePath = Path.Combine(Path.GetFullPath(OutFolder), SIM_PARAMS_FILENAME);
        using var file = new StreamWriter(filePath);
        var options = new JsonSerializerOptions { IncludeFields = true, WriteIndented = true };
        string jsonString = JsonSerializer.Serialize(simParams, options);
        file.WriteLine(jsonString);
    }
    
    public void WriteEvents(List<CNEvent> abberationList)
    {
        //TODO: Format output, talk with Tom about readable ideas
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), CN_EVENTS_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        StringBuilder abberationString = new();
        abberationString.Append(
            "Clone Name\tAbberation\tEventString\tDelta Fitness\tTotal Fitness\tNumber of mutations\n");
        foreach (var abberation in abberationList)
        {
            abberationString.Append($"{abberation.CloneName}\t" +
                                    $"{abberation.AberrationType}\t" +
                                    $"{abberation.Region}\t" +
                                    $"{Math.Round((decimal)abberation.DeltaFitness, 8).ToString(CultureInfo.InvariantCulture)}\t" +
                                    $"{Math.Round((decimal)abberation.TotalFitness, 8).ToString(CultureInfo.InvariantCulture)}\t" +
                                    $"{abberation.NrOfMutation.ToString()}\n");
        }
        outputFile.Write(abberationString.ToString());
    }
    
    public void WriteNewick(string newick)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), NEWICK_FILENAME);
        Console.WriteLine($"Writing to file {outPath}");
        using var file = new StreamWriter(outPath);
        file.Write(newick);
    }
    
    public void WriteSampleFitness(List<ProfileStats> samples)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), SAMPLE_FITNESS_FILE);
        Console.WriteLine($"Writing to file {outPath}");
        using var file = new StreamWriter(outPath);
        var myT = typeof(ProfileStats);
        var fileds = myT.GetProperties();
        var fieldNames = fileds.Select(f => f.Name).ToList();
        file.WriteLine(string.Join("\t", fieldNames));
        foreach (var sample in samples)
        {
            // Get all the field values of the record sample
            var values = fileds.Select(f => $"{f.GetValue(sample)}:F4");
            file.WriteLine(string.Join("\t", values));
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
            var newickBuild = new StringBuilder(File.ReadAllText(fileFullPath));
            return newickBuild.ToString();
        }
        catch (Exception e)
        {
            throw new Exception($"Failed to parse the file {fileFullPath}. Error {e.Message}");
        }
    }
    
    public static Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> ReadGeneLists(string folder, bool isFemale)
    {
        var geneLists = new Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>>();
        var fileMap = new Dictionary<GeneListType, string>
        {
            {GeneListType.TumorSuppressor, TSGS_TSV},
            {GeneListType.Oncogene, OGS_TSV},
            {GeneListType.Essentiality, ESSENTIALS_TSV}
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