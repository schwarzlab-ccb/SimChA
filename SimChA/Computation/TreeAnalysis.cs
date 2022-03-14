using SimChA.DataTypes;

namespace SimChA.Computation;

public class TreeAnalysis
{
    private static long SubtreeCellCount(ParentTree parentTree, TreeNode subtreeRoot)
        => subtreeRoot.Size + parentTree.Edges
            .Where(e => e.SourceId == subtreeRoot.Id)
            .Select(e => SubtreeCellCount(parentTree, parentTree.Nodes
                .Find(n => n.Id == e.TargetId)))
            .Sum();

    public static Dictionary<int, long> ComputeVAF(ParentTree parentTree)
        => parentTree.Nodes.ToDictionary(node => node.Id, node => SubtreeCellCount(parentTree, node));

    class TreeSizeData
    {
        internal int leafCount;
        internal int childCount;
    }

    private static int CountNodes(Dictionary<int, List<int>> branches, TreeSizeData data, int id, int depth)
    {
        var children = branches[id];
        if (children.Any())
        {
            data.childCount += children.Count;
            return children.Select(c => CountNodes(branches, data, c, depth + 1)).Max();
        }
        else
        {
            data.leafCount += 1;
            return depth;
        }
    }

    // Returns number of nodes, number of leafs, depth, mean child count
    public static (int, int, int, float) ComputeTreeSize(ParentTree parentTree)
    {
        TreeSizeData data = new ();
        var branches = TreeToBranches(parentTree);
        int depth = CountNodes(branches, data, parentTree.RootId, 0);
        int nodeCount = parentTree.Nodes.Count;
        int leafCount = nodeCount - data.leafCount;
        return ( nodeCount, leafCount, depth, data.childCount/(float)leafCount );
    }

    private static Dictionary<int, List<int>> TreeToBranches(ParentTree pt)
    {
        Dictionary<int, List<int>> branches = new();
        foreach (var node in pt.Nodes)
        {
            var targets = pt.Edges.Where(e => e.SourceId == node.Id).Select(e => e.TargetId).ToList();
            branches.Add(node.Id, targets);
        }
        return branches;
    }
}