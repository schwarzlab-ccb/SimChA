using SimChA.Computation;
using SimChA.DataTypes;

namespace SimChA.IO;

public class FileIO
{
    private const string DOT_FILENAME = "parent_graph.dot";
    private const string SUBCLONES_FILENAME = "subclones.out";
    private const string COPYNUMBERS_FILENAME = "copynumbers.out";
    private const string BAF_FILENAME = "baf.out";
    private const string LOGR_FILENAME = "logr.out";
    private const string POPULATIONS_DF_FILENAME = "populations.csv";
    private const string ADJACENCY_DF_FILENAME = "adjacency.csv";

    private string OutFolder { get; }

    public FileIO(string outFolder)
    {
        OutFolder = outFolder;
        Directory.CreateDirectory(outFolder);
    }

    public void WriteSubClones(IEnumerable<SubClone> subClones)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), SUBCLONES_FILENAME);
        Console.WriteLine($"Writing subclones to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        
        foreach (var subClone in subClones)
        {
            outputFile.Write(subClone);
        }
    }
    
    public void WriteParentTree(ParentTree tree)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), DOT_FILENAME);
        Console.WriteLine($"Writing parent graph to file {outPath}");
        using var outputFile = new StreamWriter(outPath);

        
        outputFile.WriteLine("Digraph SimChA {");
        foreach (var node in tree.Nodes)
        {
            outputFile.WriteLine($"\t{node.Id} [label=\"{node.Size}\"];");
        }
        foreach (var edge in tree.Edges)
        {
            outputFile.WriteLine($"\t{edge.SourceId} -> {edge.TargetId} [label=\"{edge.Distance}\"];");
        }
        outputFile.WriteLine("}");
    }

    public void WriteCopyNumbers(IEnumerable<SubClone> subClones)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), COPYNUMBERS_FILENAME);
        Console.WriteLine($"Writing CopyNumbers to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        
        outputFile.Write("sample_id\tchrom\tstart\tend\tcn_a\tcn_b\n");
        foreach (var subClone in subClones)
        {
            var copynumbers = CopyNumbers.CalcCopyNumbers(subClone.Karyotype);
            outputFile.Write(CopyNumbers.ToTSV(copynumbers, subClone.CloneId.ToString(), false) + "\n");
        }
    }

    public void WriteMullerDataFrames(IEnumerable<SubClone> subClones)
    {   
        string popPath = Path.Combine(Path.GetFullPath(OutFolder), POPULATIONS_DF_FILENAME);
        string adjPath = Path.Combine(Path.GetFullPath(OutFolder), ADJACENCY_DF_FILENAME);
        Console.WriteLine($"Writing population to file {popPath}, adjacency to file {adjPath}");
        
        using var popFile = new StreamWriter(popPath);
        popFile.WriteLine("Generation,Identity,Population");
        foreach (var subClone in subClones)
        {
            for (int i = 0; i < subClone.Generations.Count; i++)
            {
                popFile.WriteLine($"{subClone.FirstGen + i},{subClone.CloneId},{subClone.Generations[i]}");
            }
        }
        
        using var adjFile = new StreamWriter(adjPath);
        adjFile.WriteLine("Parent,Identity");
        foreach (var subClone in subClones)
        {
            adjFile.WriteLine($"{subClone.ParentId},{subClone.CloneId}");
        }
    }

    public void WriteRawData(IEnumerable<SubClone> subClones, List<SNP> snps, bool isFemale)
    {
        var outputbaf = new List<string>();
        outputbaf.Add("\tchrom\tpos");

        foreach (var snp in snps)
        {
            outputbaf.Add($"{snp.Id}\t{snp.Chrom}\t{snp.Pos}");
        }
        var outputlogr = outputbaf.Select(x => x.Clone()).ToList();
        foreach (var subClone in subClones)
        {
            var copyNumbers = CopyNumbers.CalcCopyNumbers(subClone.Karyotype);
            var rawdata = RawData.CalcSingleSubclone(copyNumbers, snps, isFemale);
            outputbaf[0] += $"\t{subClone.CloneId}";
            outputlogr[0] += $"\t{subClone.CloneId}";
            for (int i = 0; i < rawdata.Count; i++)
            {
                outputbaf[i + 1] += $"\t{rawdata[i].Baf}";
                outputlogr[i + 1] += $"\t{rawdata[i].LogR}";
            }
        }

        string outPathBAF = Path.Combine(Path.GetFullPath(OutFolder), BAF_FILENAME);
        Console.WriteLine($"Writing BAF to file {outPathBAF}");
        using var outputFileBAF = new StreamWriter(outPathBAF);
        outputFileBAF.Write(string.Join("\n", outputbaf) + "\n");

        string outPathLogR = Path.Combine(Path.GetFullPath(OutFolder), LOGR_FILENAME);
        Console.WriteLine($"Writing LogR to file {outPathLogR}");
        using var outputFileLogR = new StreamWriter(outPathLogR);
        outputFileLogR.Write(string.Join("\n", outputlogr) + "\n");
    }
    public void WriteRawData(List<SNPData> rawData, int subcloneId)
    {
        string outPathBAF = Path.Combine(Path.GetFullPath(OutFolder), $"{subcloneId}_{BAF_FILENAME}");
        Console.WriteLine($"Writing BAF to file {outPathBAF}");
        using var outputFileBAF = new StreamWriter(outPathBAF);
        outputFileBAF.Write(RawData.PrintBAF(rawData) + "\n");

        string outPathLogR = Path.Combine(Path.GetFullPath(OutFolder), $"{subcloneId}_{LOGR_FILENAME}");
        Console.WriteLine($"Writing LogR to file {outPathLogR}");
        using var outputFileLogR = new StreamWriter(outPathLogR);
        outputFileLogR.Write(RawData.PrintLogR(rawData) + "\n");
    }
}