using System.Globalization;
using System.Text.Json;
using SimChA.Computation;
using SimChA.DataTypes;

namespace SimChA.IO;

public class FileIO
{
    private const string DOT_FILENAME = "parent_graph.dot";
    private const string SUBCLONES_FILENAME = "subclones.out";
    // private const string COPYNUMBERS_FILENAME = "copynumbers.out";
    private const string BAF_FILENAME = "baf.out";
    private const string LOGR_FILENAME = "logr.out";
    private const string POPULATIONS_DF_FILENAME = "populations.csv";
    private const string ADJACENCY_DF_FILENAME = "parent_tree.csv";
    private const string SIM_PARAMS_FILENAME = "sim_params.json";
    private const string CCF_FILENAME = "ccf.csv";
    private const string SUMMARY_FILENAME = "summary.csv";
    private string Timestamp { get; }
    private string RootFolder { get; }
    private string ExperimentFolder { get; }

    public FileIO(string rootFolder)
    {
        Timestamp =  DateTime.Now.ToString("yy_MM_dd_hh_mm_ss");
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

        ExperimentFolder = Path.Join(RootFolder, Timestamp);
        Directory.CreateDirectory(ExperimentFolder);
    }

    public void WriteSubClones(IEnumerable<SubClone> subClones)
    {
        string outPath = Path.Combine(Path.GetFullPath(RootFolder), SUBCLONES_FILENAME);
        using var outputFile = new StreamWriter(outPath);
        
        foreach (var subClone in subClones)
        {
            outputFile.Write(subClone);
        }
    }
    
    public void WriteParentTree(ParentTree tree)
    {
        string outPath = Path.Combine(Path.GetFullPath(RootFolder), DOT_FILENAME);
        using var outputFile = new StreamWriter(outPath);
        
        outputFile.WriteLine("Digraph SimChA {");
        foreach (var node in tree.Nodes)
        {
            double size = Math.Round(.25 * (1 + Math.Log(1 + node.Size)), 2);
            outputFile.WriteLine($"\t{node.Id} [label=\"{node.Id}:{node.Size}\", width={size}, height={size*.6}];");
        }
        foreach (var edge in tree.Edges)
        {
            outputFile.WriteLine($"\t{edge.SourceId} -> {edge.TargetId} [label=\"{edge.Distance}\"];");
        }
        outputFile.WriteLine("}");
    }

    // public void WriteCopyNumbers(IEnumerable<SubClone> subClones)
    // {
    //     string outPath = Path.Combine(Path.GetFullPath(OutFolder), COPYNUMBERS_FILENAME);
    //     Console.WriteLine($"Writing CopyNumbers to file {outPath}");
    //     using var outputFile = new StreamWriter(outPath);
    //     
    //     outputFile.Write("sample_id\tchrom\tstart\tend\tcn_a\tcn_b\n");
    //     foreach (var subClone in subClones)
    //     {
    //         var copynumbers = CopyNumbers.CalcCopyNumbers(subClone.Karyotype);
    //         outputFile.Write(CopyNumbers.ToTSV(copynumbers, subClone.CloneId.ToString(), false) + "\n");
    //     }
    // }

    public void WriteMullerDataFrames(IEnumerable<SubClone> subClones, ParentTree tree)
    {   
        string popPath = Path.Combine(Path.GetFullPath(RootFolder), POPULATIONS_DF_FILENAME);
        string adjPath = Path.Combine(Path.GetFullPath(RootFolder), ADJACENCY_DF_FILENAME);

        using var popFile = new StreamWriter(popPath);
        popFile.WriteLine("Id,Step,Pop");
        foreach (var subClone in subClones)
        {
            int start = subClone.FirstGen;
            int end = subClone.LastGen;
            for (int gen = start; gen < end ; gen++)
            {
                int alive = subClone.AliveAtGen(gen);
                if (alive > 0)
                {
                    popFile.WriteLine($"{subClone.CloneId},{gen},{alive}");
                }
            }
        }

        using var adjFile = new StreamWriter(adjPath);
        adjFile.WriteLine("ParentId,ChildId");

        foreach (var edge in tree.Edges)
        {
            adjFile.WriteLine($"\t{edge.SourceId},{edge.TargetId}");
        }
    }

    public void WriteCCF(Dictionary<int, long> vaf, long totalSize)
    {
        string outPath = Path.Combine(Path.GetFullPath(RootFolder), CCF_FILENAME);
        using var outputFile = new StreamWriter(outPath);
        outputFile.WriteLine("id,pop,ccf");
        foreach ((int id, long pop) in vaf)
        {
            outputFile.WriteLine($"{id},{pop},{(float) pop / totalSize}");
        }
    }
    
    // public void WriteRawData(Random rnd, IEnumerable<SubClone> subClones, List<SNP> snps, bool isFemale)
    // {
    //     var outputbaf = new List<string>();
    //     outputbaf.Add("\tchrom\tpos");
    //
    //     foreach (var snp in snps)
    //     {
    //         outputbaf.Add($"{snp.Id}\t{snp.Chrom}\t{snp.Pos}");
    //     }
    //     var outputlogr = outputbaf.Select(x => x.Clone()).ToList();
    //     foreach (var subClone in subClones)
    //     {
    //         var copyNumbers = CopyNumbers.CalcCopyNumbers(subClone.Karyotype);
    //         var rawdata = SNPMetrics.CalcSingleSubClone(rnd, copyNumbers, snps, isFemale);
    //         outputbaf[0] += $"\t{subClone.CloneId}";
    //         outputlogr[0] += $"\t{subClone.CloneId}";
    //         for (int i = 0; i < rawdata.Count; i++)
    //         {
    //             outputbaf[i + 1] += $"\t{rawdata[i].Baf}";
    //             outputlogr[i + 1] += $"\t{rawdata[i].LogR}";
    //         }
    //     }
    //
    //     string outPathBAF = Path.Combine(Path.GetFullPath(OutFolder), BAF_FILENAME);
    //     Console.WriteLine($"Writing BAF to file {outPathBAF}");
    //     using var outputFileBAF = new StreamWriter(outPathBAF);
    //     outputFileBAF.Write(string.Join("\n", outputbaf) + "\n");
    //
    //     string outPathLogR = Path.Combine(Path.GetFullPath(OutFolder), LOGR_FILENAME);
    //     Console.WriteLine($"Writing LogR to file {outPathLogR}");
    //     using var outputFileLogR = new StreamWriter(outPathLogR);
    //     outputFileLogR.Write(string.Join("\n", outputlogr) + "\n");
    // }
    
    public void WriteRawData(List<SNPData> rawData, int subcloneId)
    {
        string outPathBAF = Path.Combine(Path.GetFullPath(RootFolder), $"{subcloneId}_{BAF_FILENAME}");
        using var outputFileBAF = new StreamWriter(outPathBAF);
        outputFileBAF.Write(SNPMetrics.PrintBAF(rawData) + "\n");

        string outPathLogR = Path.Combine(Path.GetFullPath(RootFolder), $"{subcloneId}_{LOGR_FILENAME}");
        using var outputFileLogR = new StreamWriter(outPathLogR);
        outputFileLogR.Write(SNPMetrics.PrintLogR(rawData) + "\n");
    }

    public void WriteSimParams(SimParams simParams)
    {
        string filePath = Path.Combine(Path.GetFullPath(ExperimentFolder), SIM_PARAMS_FILENAME);
        using var file = new StreamWriter(filePath);
        var options = new JsonSerializerOptions { IncludeFields = true, WriteIndented = true };
        string jsonString = JsonSerializer.Serialize(simParams, options);
        file.WriteLine(jsonString);
    }

    public void CreateSummary()
    {
        string filePath = Path.Combine(Path.GetFullPath(ExperimentFolder), SUMMARY_FILENAME);
        using var file = new StreamWriter(filePath);
        file.WriteLine(ResultSummary.Header());
    }

    public void AddToSummary(ResultSummary resultSummary)
    {
        string filePath = Path.Combine(Path.GetFullPath(ExperimentFolder), SUMMARY_FILENAME);
        using var file = new StreamWriter(filePath, true);
        file.WriteLine(resultSummary.ToString());
    }

    public void StoreCopy(int runId)
    {
        string copyFolder = Path.Join(ExperimentFolder, runId.ToString());
        Directory.CreateDirectory(copyFolder);
        
        foreach (var file in new DirectoryInfo(RootFolder).GetFiles())
        {
            file.CopyTo(Path.Join(copyFolder, file.Name));
        }
    }

    public void CopySummary()
    {
        File.Copy(
            Path.Combine(Path.GetFullPath(ExperimentFolder), SUMMARY_FILENAME), 
            Path.Combine(Path.GetFullPath(RootFolder), SUMMARY_FILENAME));
        
        File.Copy(
            Path.Combine(Path.GetFullPath(ExperimentFolder), SIM_PARAMS_FILENAME), 
            Path.Combine(Path.GetFullPath(RootFolder), SIM_PARAMS_FILENAME));
    }
}