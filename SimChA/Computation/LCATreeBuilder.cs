using SimChA.DataTypes;

namespace SimChA.Computation;

public static class LcaTreeBuilder
{
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
    public static ParentTree BuildTree(List<Clone> selection)
    {
        List<TreeNode> nodes = new();
        List<TreeEdge> edges = new();

        foreach (var subClone in selection)
        {
            nodes.Add(new TreeNode {Id = subClone.CloneId, Name = subClone.Name});
            foreach(var childId in subClone.ChildrenIDs){
                edges.Add(new TreeEdge{
                    Distance = selection[childId].MutCount, SourceId = subClone.Name, TargetId = selection[childId].Name});
            }
        }

        return new ParentTree {RootName = selection[0].Name, Nodes = nodes, Edges = edges};
    }
}