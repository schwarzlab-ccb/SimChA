using SimChA.DataTypes;

namespace SimChA.IO;

public class DotGraph
{
    public static void WriteGraph(string path, List<SubClone> subClones)
    {
        Console.WriteLine($"Writing DoT graph to file {path}");
        using var outputFile = new StreamWriter(path);
        outputFile.WriteLine("Digraph SimChA {");
        foreach (var subClone in subClones.Where(sc => sc.AliveCount > 10))
        {
            outputFile.WriteLine($"\t{subClone.CloneId} [label=\"{subClone.AliveCount}\"];");
            if (subClone.ParentId >= 0)
            {
                outputFile.WriteLine($"\t{subClone.ParentId} -> {subClone.CloneId};");
            }
        }
        outputFile.WriteLine("}");
    }
}