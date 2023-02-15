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
    private const string TSV_FILENAME = "abberations.tsv";
    private const string ESSENTIALS_TSV = "essentials.tsv";
    private const string OGS_TSV = "ogs.tsv";
    private const string TSGS_TSV = "tsgs.tsv";
    private const string SAMPLE_FITNESS_FILE = "sample_fitness.tsv";

    public FileIO(string rootFolder)
    {
        Timestamp = DateTime.Now.ToString("yy_MM_dd_HH_mm_ss");
        CultureInfo.DefaultThreadCurrentCulture = CultureInfo.InvariantCulture;

        RootFolder = rootFolder;
        if (Directory.Exists(RootFolder))
        {
            foreach (var file in new DirectoryInfo(RootFolder).GetFiles())
            {
                file.Delete();
            }
        }
        else
        {
            Directory.CreateDirectory(RootFolder);
        }

        if (IsRepeated)
        {
            ExperimentFolder = Path.Join(RootFolder, Timestamp);
            Directory.CreateDirectory(ExperimentFolder);
        }
        else
        {
            ExperimentFolder = RootFolder;
        }
    }

    private string Timestamp { get; }
    private string RootFolder { get; }
    private string ExperimentFolder { get; }
    private bool IsRepeated { get; }

    public void WriteClones(IEnumerable<Clone> subClones)
    {
        string outPath = Path.Combine(Path.GetFullPath(RootFolder), SUBCLONES_FILENAME);
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
        string outPath = Path.Combine(Path.GetFullPath(RootFolder), COPYNUMBERS_FILENAME);
        Console.WriteLine($"Writing CopyNumbers to file {outPath}");
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
        string filePath = Path.Combine(Path.GetFullPath(ExperimentFolder), SIM_PARAMS_FILENAME);
        using var file = new StreamWriter(filePath);
        var options = new JsonSerializerOptions { IncludeFields = true, WriteIndented = true };
        string jsonString = JsonSerializer.Serialize(simParams, options);
        file.WriteLine(jsonString);
    }

    public static SimParams SimParamsFromFile(string filePath)
    {
        string fileFullPath = Path.GetFullPath(filePath);
        if (!File.Exists(fileFullPath))
        {
            throw new Exception($"Configuration file {fileFullPath} does not exist");
        }

        try
        {
            string serializedJSON = File.ReadAllText(fileFullPath);
            var options = new JsonSerializerOptions { IncludeFields = true };
            var simParams = JsonSerializer.Deserialize<SimParams>(serializedJSON, options);
            if (simParams.Seed < 0)
            {
                return simParams with { Seed = new Random().Next() };
            }
            return simParams;
        }
        catch (Exception e)
        {
            throw new Exception($"Failed to read simulation params from the file {fileFullPath}. Error {e.Message}");
        }
    }

    public static string GetStringFromNewick(string newickFile)
    {
        string fileFullPath = Path.GetFullPath(newickFile);
        if (!File.Exists(fileFullPath))
        {
            throw new Exception($"Newick file {fileFullPath} does not exist");
        }
        try
        {
            var newickBuild = new StringBuilder(File.ReadAllText(fileFullPath));
            return newickBuild.ToString();
        }
        catch (Exception e)
        {
            throw new Exception($"Failed to read newick file from the file {fileFullPath}. Error {e.Message}");
        }
    }

    public void WriteTSV(List<Abberation> abberationList)
    {
        //TODO: Format output, talk with Tom about readable ideas
        string filePath = Path.Combine(Path.GetFullPath(RootFolder), TSV_FILENAME);
        using var outputFile = new StreamWriter(filePath);
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
            if (!File.Exists(filePath))
            {
                throw new Exception($"Required file {filename} not found in {folder} directory.");
            }
            string fileContent = File.ReadAllText(Path.GetFullPath(filePath));
            geneLists[key] = ReadGenesFromFile(fileContent, isFemale);
        }
        return geneLists;
    }

    public static Dictionary<ChrNo, List<Gene>> ReadGenesFromFile(string fileContent, bool isFemale)
    {
        // Pre-initialization
        var noEnum = Enum.GetValues(typeof(ChrNo)).Cast<ChrNo>().ToList();
        var geneList = noEnum.ToDictionary(c => c, c => new List<Gene>());
        string[] genesFromFile = fileContent.Split('\n', StringSplitOptions.RemoveEmptyEntries);
        foreach (string geneFromFile in genesFromFile)
        {
            string[] genString = geneFromFile.Split('\t');
            //Don't include Y chromosome in genes list if clone is female
            if (isFemale && (ChrNo)Enum.Parse(typeof(ChrNo), genString[2]) == ChrNo.chrY)
            {
                continue;
            }
            string name = genString[3];
            double fitness = double.Parse(genString[4], CultureInfo.InvariantCulture.NumberFormat);
            var chrNum = (ChrNo)Enum.Parse(typeof(ChrNo), genString[0]);
            var chrID = new ChrID(chrNum, isFemale);
            // Convert to zero-based [start, end) index 
            var region = new Region(int.Parse(genString[1]) - 1, int.Parse(genString[2]), chrID);
            var gene = new Gene(name, region, fitness);
            geneList[chrNum].Add(gene);
        }
        return geneList;
    }

    public void WriteSampleFitness(Dictionary<string, double> fitness)
    {
        string filePath = Path.Combine(Path.GetFullPath(ExperimentFolder), SAMPLE_FITNESS_FILE);
        using var file = new StreamWriter(filePath);
        file.WriteLine("Sample\tFitness");
        foreach (var (sample, fit) in fitness)
        {
            file.WriteLine($"{sample}\t{fit}");
        }
    }

    public static Dictionary<string, Karyotype> ReadCopyNumbers(string cnaProfile)
    {
        string fileFullPath = Path.GetFullPath(cnaProfile);
        Dictionary<string, Karyotype> result = new();
        if (!File.Exists(fileFullPath))
        {
            throw new Exception($"CNA profile file {fileFullPath} does not exist");
        }
        try
        {
            var cnaFile = new StreamReader(fileFullPath);
            string lastSample = "";
            var regionsA = new List<Region>();
            var regionsB = new List<Region>();
            cnaFile.ReadLine(); // Skip header
            while (cnaFile.ReadLine() is { } line)
            {
                string[] lineSplit = line.Split('\t');
                string sample = lineSplit[0];
                if (sample != lastSample)
                {
                    if (regionsA.Any() || regionsB.Any())
                    {
                        var haplotypes = new List<Contig> {new(regionsA), new(regionsB)};
                        result[sample] = new Karyotype(haplotypes);
                        regionsA.Clear();
                        regionsB.Clear();
                    }
                    lastSample = sample;
                }
                var chrNum = (ChrNo)Enum.Parse(typeof(ChrNo), lineSplit[1]);
                int start = int.Parse(lineSplit[2]) - 1;
                int end = int.Parse(lineSplit[3]);
                int cnA = int.Parse(lineSplit[4]);
                int cnB = int.Parse(lineSplit[5]);
                for (int i = 0; i < cnA; i++)
                {
                    regionsA.Add(new Region(start, end, new ChrID(chrNum, true)));
                }
                for (int i = 0; i < cnB; i++)
                {
                    regionsB.Add(new Region(start, end, new ChrID(chrNum, false)));
                }
            }
        }
        catch (Exception e)
        {
            throw new Exception($"Failed to read newick file from the file {fileFullPath}. Error {e.Message}");
        }
        return result;
    }
}