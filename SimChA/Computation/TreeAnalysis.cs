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
}