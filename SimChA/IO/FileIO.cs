using System.Globalization;
using System.Text;
using System.Text.Json;
using SimChA.Computation;
using SimChA.DataTypes;

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

        StringBuilder clonesString = new StringBuilder();
        foreach (var subClone in subClones)
        {
            clonesString.AppendLine(subClone.ToString());
        }

        outputFile.Write(clonesString.ToString());
    }

    public void WriteParentTree(ParentTree tree)
    {
        string outPath = Path.Combine(Path.GetFullPath(RootFolder), DOT_FILENAME);
        using var outputFile = new StreamWriter(outPath);

        StringBuilder parentTreeString = new StringBuilder();
        parentTreeString.AppendLine("Digraph SimChA {");
        foreach (var node in tree.Nodes)
        {
            double size = Math.Round(.25 * (1 + Math.Log(1 + 0)), 2);
            parentTreeString.AppendLine($"\t{node.Id} [label=\"{node.Name}\", width={size}, height={size * .6}];");
        }

        foreach (var edge in tree.Edges)
        {
            parentTreeString.AppendLine($"\t{edge.SourceId} -> {edge.TargetId} [label=\"{edge.Distance}\"];");
        }

        parentTreeString.AppendLine("}");
        outputFile.Write(parentTreeString.ToString());
    }

    public void WriteNewickFile(List<Clone> clones)
    {
        string outPath = Path.Combine(Path.GetFullPath(RootFolder), NEWICK_FILENAME);
        using var outputFile = new StreamWriter(outPath);
        StringBuilder newickString = new StringBuilder(";");
        newickString.Insert(0, IterateClones(clones[0], clones));
        outputFile.WriteLine(newickString);
    }

    private static string IterateClones(Clone clone, List<Clone> clones)
    {
        var newickString = new StringBuilder(":" + clone.DistToParent);
        newickString.Insert(0, clone.Name);
        newickString.Insert(0, clone.ChildrenIDs.Count > 0 ? ")":"");
        foreach(int cloneId in clone.ChildrenIDs){
            newickString.Insert(0, IterateClones(clones[cloneId], clones));
            newickString.Insert(0, cloneId != clone.ChildrenIDs.Last() ? "," : "");
        }
        newickString.Insert(0, clone.ChildrenIDs.Count > 0 ? "(":"");
        return newickString.ToString();
    }

    public void WriteCopyNumbers(IEnumerable<Clone> subClones)
    {
        string outPath = Path.Combine(Path.GetFullPath(RootFolder), COPYNUMBERS_FILENAME);
        Console.WriteLine($"Writing CopyNumbers to file {outPath}");
        using var outputFile = new StreamWriter(outPath);

        StringBuilder copyNumbersString = new StringBuilder();
        copyNumbersString.Append("sample_id\tchrom\tstart\tend\tcn_a\tcn_b\n");

        foreach (var subClone in subClones)
        {
            var copynumbers = CopyNumbers.CalcCopyNumbers(subClone.Karyotype);
            copyNumbersString.Append((CopyNumbers.ToTSV(copynumbers, subClone.Name.ToString(), false) + "\n"));
        }
        outputFile.Write(copyNumbersString.ToString());
    }
    
    public void WriteRawData(Random rnd, IEnumerable<Clone> subClones, List<SNP> snps, bool isFemale)
    {
        var outputbaf = new List<string>();
        outputbaf.Add("\tchrom\tpos");

        outputbaf.AddRange(snps.Select(snp => $"{snp.Id}\t{snp.Chrom}\t{snp.Pos}"));

        var outputlogr = outputbaf.Select(x => x.Clone()).ToList();
        foreach (var subClone in subClones)
        {
            var copyNumbers = CopyNumbers.CalcCopyNumbers(subClone.Karyotype);
            var rawdata = SNPMetrics.CalcSingleSubClone(rnd, copyNumbers, snps, isFemale);
            outputbaf[0] += $"\t{subClone.CloneId}";
            outputlogr[0] += $"\t{subClone.CloneId}";
            for (int i = 0; i < rawdata.Count; i++)
            {
                outputbaf[i + 1] += $"\t{rawdata[i].Baf}";
                outputlogr[i + 1] += $"\t{rawdata[i].LogR}";
            }
        }

        string outPathBAF = Path.Combine(Path.GetFullPath(RootFolder), BAF_FILENAME);
        Console.WriteLine($"Writing BAF to file {outPathBAF}");
        using var outputFileBAF = new StreamWriter(outPathBAF);
        outputFileBAF.Write(string.Join("\n", outputbaf) + "\n");

        string outPathLogR = Path.Combine(Path.GetFullPath(RootFolder), LOGR_FILENAME);
        Console.WriteLine($"Writing LogR to file {outPathLogR}");
        using var outputFileLogR = new StreamWriter(outPathLogR);
        outputFileLogR.Write(string.Join("\n", outputlogr) + "\n");
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
        StringBuilder abberationString = new StringBuilder();
        abberationString.Append($"Clone Name\tAbberation\tEventString\tdelta Fitness\tTotal Fitness\tNumber of mutations\n");
        foreach(Abberation abberation in abberationList){
            abberationString.Append($"{abberation.CloneName}\t" + 
                                    $"{abberation.AbberationEnum}\t" +
                                    $"{abberation.Region}\t" +
                                    $"{Math.Round((decimal)abberation.DeltaFitness,8).ToString()}\t" +
                                    $"{Math.Round((decimal)abberation.TotalFitness,8).ToString()}\t" +
                                    $"{abberation.NrOfMutation.ToString()}\n");
        }
        outputFile.Write(abberationString.ToString());
    }

    public static (Dictionary<ChromNum, List<Gene>>, Dictionary<ChromNum, List<Gene>>) ReadGenes(string folder, bool isFemale)
    {
        Dictionary<ChromNum, List<Gene>> tsgOgList = new();
        Dictionary<ChromNum, List<Gene>> essentialList = new();
        
        Directory.GetFiles(Path.GetFullPath("./"));
        if(!File.Exists(Path.Combine(folder, TSGS_TSV)) 
           || !File.Exists(Path.Combine(folder, OGS_TSV)) 
           || !File.Exists(Path.Combine(folder, ESSENTIALS_TSV)))
        {
            throw new Exception($"Required files not found in {folder} directory.");
        }
        ReadGenesFromFile(Path.Combine(folder, TSGS_TSV), true, tsgOgList, isFemale);
        ReadGenesFromFile(Path.Combine(folder, OGS_TSV), false, tsgOgList, isFemale);
        ReadGenesFromFile(Path.Combine(folder, ESSENTIALS_TSV), true, essentialList, isFemale);        
        return (tsgOgList, essentialList);
    }

    private static void ReadGenesFromFile(string file, bool negative, Dictionary<ChromNum, List<Gene>> genes, bool isFemale)
    {
        string fileContent = File.ReadAllText(Path.GetFullPath(file));
        string[] genesFromFile = fileContent.Split('\n');
        foreach(string geneFromFile in genesFromFile)
        {   
            if(geneFromFile != "")
            {
                string[] genString = geneFromFile.Split('\t');
                //Don't include Y chromosome in genes list if clone is female
                if(isFemale && (ChromNum)System.Enum.Parse(typeof(ChromNum), genString[2]) == ChromNum.chrY)
                {
                    continue;
                }
                string name = genString[0];
                float fitness = float.Parse(genString[1]);
                var chromNum = (ChromNum) Enum.Parse(typeof(ChromNum), genString[2]);
                var chromID = new ChromID(chromNum, false);
                string[] splits = genString[4].Split('\r');
                var region = new Region(int.Parse(genString[3]), int.Parse(splits[0]), chromID);
                var gene = new Gene(name, region, negative ? -fitness : fitness);
                if(genes.ContainsKey(chromNum))
                {
                    genes[chromNum].Add(gene);
                }
                else
                {
                    genes.Add(chromNum, new List<Gene> {gene});
                }
            }

        }
    }
}