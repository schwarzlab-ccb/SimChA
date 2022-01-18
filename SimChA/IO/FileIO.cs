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

    private string OutFolder { get; }

    public FileIO(string outFolder)
    {
        OutFolder = outFolder;
        Directory.CreateDirectory(outFolder);
    }
    
    public void WriteSubClones(IEnumerable<SubClone> subClones, int cutOff)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), SUBCLONES_FILENAME);
        Console.WriteLine($"Writing subclones to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        
        foreach (var subClone in subClones.Where(subClone => subClone.AliveCount >= cutOff))
        {
            outputFile.Write(subClone);
        }
    }
    
    public void WriteParentGraph(IEnumerable<SubClone> subClones, int cutOff)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), DOT_FILENAME);
        Console.WriteLine($"Writing parent graph to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        
        outputFile.WriteLine("Digraph SimChA {");
        foreach (var subClone in subClones.Where(sc => sc.AliveCount > cutOff))
        {
            outputFile.WriteLine($"\t{subClone.CloneId} [label=\"{subClone.AliveCount}\"];");
            if (subClone.ParentId >= 0)
            {
                outputFile.WriteLine($"\t{subClone.ParentId} -> {subClone.CloneId};");
            }
        }
        outputFile.WriteLine("}");
    }

    public void WriteCopyNumbers(IEnumerable<SubClone> subClones, int cutOff)
    {
        string outPath = Path.Combine(Path.GetFullPath(OutFolder), COPYNUMBERS_FILENAME);
        Console.WriteLine($"Writing CopyNumbers to file {outPath}");
        using var outputFile = new StreamWriter(outPath);
        
        outputFile.Write("sample_id\tchrom\tstart\tend\tcn_a\tcn_b\n");
        foreach (var subClone in subClones.Where(sc => sc.AliveCount > cutOff))
        {
            var copynumbers = CopyNumbers.CalcCopyNumbers(subClone.Karyotype);
            outputFile.Write(CopyNumbers.ToTSV(copynumbers, subClone.CloneId.ToString(), false) + "\n");
        }
    }

    public void WriteRawData(List<SNPData> rawData)
    {
        string outPathBAF = Path.Combine(Path.GetFullPath(OutFolder), BAF_FILENAME);
        Console.WriteLine($"Writing BAF to file {outPathBAF}");
        using var outputFileBAF = new StreamWriter(outPathBAF);
        outputFileBAF.Write(RawData.PrintBAF(rawData) + "\n");

        string outPathLogR = Path.Combine(Path.GetFullPath(OutFolder), LOGR_FILENAME);
        Console.WriteLine($"Writing LogR to file {outPathLogR}");
        using var outputFileLogR = new StreamWriter(outPathLogR);
        outputFileLogR.Write(RawData.PrintLogR(rawData) + "\n");
    }
}