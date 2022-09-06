using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.IO;

public static class Newick
{
    private static Clone CreateNodes(string newickNode, int parentId, bool isFemale)
    {
        string[] cloneString = newickNode.Split(':');
        // TODO: split below in individual assignments
        var clone = new Clone(int.Parse(cloneString[0].Split('-')[0]), parentId, int.Parse(cloneString[1]),
            int.Parse(cloneString[0].Split('-')[1]), new Karyotype(isFemale));
        return clone;
    }

    public static List<Clone> ParseNewickTree(string[] newickString, bool isFemale)
    {
        List<Clone> clones = new();
        var parentIds = new List<int> { -1 };
        bool rootSet = false;
        // TODO: Multiple code repetitions below, fix
        for (int i = 0; i < newickString.Length; i++)
        {
            switch (newickString[i])
            {
                case "(":
                    if (newickString[i - 1] == "")
                    {
                        parentIds = parentIds.Where(p => p != parentIds.Last()).ToList();
                        break;
                    }

                    clones.Add(CreateNodes(newickString[i - 1], parentIds.Last(), isFemale));
                    parentIds = parentIds.Where(p => p != parentIds.Last()).ToList();
                    break;
                case ")":
                    if (rootSet)
                    {
                        clones.Add(CreateNodes(newickString[i - 1], parentIds.Last(), isFemale));
                        parentIds.Add(int.Parse(newickString[i - 1].Split('-')[0]));
                    }

                    break;
                case ",":
                    if (!rootSet)
                    {
                        clones.Add(CreateNodes(newickString[i - 1], parentIds.Last(), isFemale));
                        parentIds.Add(int.Parse(newickString[i - 1].Split('-')[0]));
                        rootSet = true;
                    }
                    else
                    {
                        clones.Add(CreateNodes(newickString[i - 1], parentIds.Last(), isFemale));
                    }

                    break;
            }
        }
        return clones;
    }
}