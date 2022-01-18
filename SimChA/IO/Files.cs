using SimChA.DataTypes;

namespace SimChA.IO;

public class Files
{
    private const string DOT_FILENAME = "parent_graph.dot";
    private const string SUBCLONES_FILENAME = "subclones.txt";
    private const string COPYNUMBERS_FILENAME = "copynumbers.out";

    private string OutFolder { get; }

    public Files(string outFolder)
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
            var copynumbers = new CopyNumbers(subClone.Karyotype);
            outputFile.Write(copynumbers.ToTSV(subClone.CloneId.ToString(), false));
        }
    }
}