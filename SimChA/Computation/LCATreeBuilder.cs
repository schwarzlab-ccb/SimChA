using SimChA.DataTypes;

namespace SimChA.Computation;

public static class LCATreeBuilder
{
    public static bool isNewick { get; set; } = false;

    private static TreeEdge FindEdge(Dictionary<int, int> parentMap, List<Clone> selection, List<int> internalNodes,
        int id)
    {
        int dist = 0;
        int source = id;
        do
        {
            dist++;
            source = parentMap[source];
        } while (source != -1 && selection.All(sc => sc.CloneId != source) && internalNodes.All(n => n != source));

        return new TreeEdge { Distance = dist, SourceId = source, TargetId = id };
    }

    private static List<int> FindInternalNodes(Dictionary<int, int> parentMap, List<Clone> selection)
    {
        Dictionary<int, int> internalNodes = new();

        foreach (var subClone in selection)
        {
            int curNode = parentMap[subClone.CloneId];
            while (selection.All(sc => sc.CloneId != curNode) && curNode != -1)
            {
                if (internalNodes.ContainsKey(curNode))
                {
                    internalNodes[curNode]++;
                    break;
                }

                internalNodes[curNode] = 0;
                curNode = parentMap[curNode];
            }
        }

        return internalNodes.Where(n => n.Value > 0 || n.Key == 0).Select(n => n.Key).ToList();
    }

    // Construct a parent tree with lowest common ancestor (LCA) for each pair of children
    public static ParentTree Builtree(IEnumerable<Clone> allSubClones, List<Clone> selection)
    {
        var parentMap = ConnectedTreeBuilder.CreateParentMap(allSubClones);
        var internalNodes = FindInternalNodes(parentMap, selection);

        List<TreeNode> nodes = new();
        List<TreeEdge> edges = new();

        foreach (var subClone in selection)
        {
            nodes.Add(new TreeNode { Id = subClone.CloneId, Size = subClone.CellCount });
            if (!isNewick)
            {
                edges.Add(FindEdge(parentMap, selection, internalNodes, subClone.CloneId));
            }
            else
            {
                edges.Add(new TreeEdge
                    { Distance = subClone.MutCount, SourceId = subClone.ParentId, TargetId = subClone.CloneId });
            }
        }

        foreach (int internalNode in internalNodes)
        {
            nodes.Add(new TreeNode { Id = internalNode, Size = 0 });
            if (!isNewick)
            {
                edges.Add(FindEdge(parentMap, selection, internalNodes, internalNode));
            }
            else
            {
                var currentPos = allSubClones.Where(id => id.CloneId == internalNode).First();
                edges.Add(new TreeEdge
                    { Distance = currentPos.MutCount, SourceId = currentPos.ParentId, TargetId = internalNode });
            }
        }

        return new ParentTree { RootId = 0, Nodes = nodes, Edges = edges.Where(e => e.TargetId != 0).ToList() };
    }
}