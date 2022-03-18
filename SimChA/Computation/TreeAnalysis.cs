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
        int leafCount = data.leafCount;
        float branching = data.childCount / (float)leafCount;
        return ( nodeCount, leafCount, depth, branching );
    }

    public static float ComputeTreeBalance(ParentTree parentTree)
    {
        float treeBalance = 0;
        long Sdash_i_sum = 0;
        var subtreeCount = ComputeVAF(parentTree);
        var branches = TreeToBranches(parentTree);

        foreach (var node in parentTree.Nodes.Where(n => branches[n.Id].Count() >= 2))
        {
            int nChildren = branches[node.Id].Count();
            long S_i = subtreeCount[node.Id];
            long Sdash_i = S_i - node.Size;
            Sdash_i_sum += Sdash_i;

            float W_i = branches[node.Id].Select(b => (float)subtreeCount[b] / Sdash_i)
                                         .Where(p => p > 0)
                                         .Select(p => -1 * p * (float)Math.Log(p) / (float)Math.Log(nChildren))
                                         .Sum();

            treeBalance += Sdash_i * Sdash_i / S_i * W_i;
        }

        return treeBalance / Sdash_i_sum;
    }

    public static float ComputeClonalDiversity(ParentTree parentTree)
    {
        long totalPop = SubtreeCellCount(parentTree, parentTree.Nodes.Where(n => n.Id == parentTree.RootId).First());
        float clonalDiversity = 1 / parentTree.Nodes.Where(node => node.Size > 0)
                                                .Select(node => (float)Math.Pow(((float)node.Size / totalPop), 2))
                                                .Sum();

        return clonalDiversity;
    }

    public static float ComputeClonalDiversityFiltered(List<SubClone> subClones)
    {
        long totalPop = subClones.Select(clone => clone.AliveCount).Sum();
        float clonalDiversity = 1 / subClones.Select(clone => (float)Math.Pow(((float)clone.AliveCount / totalPop), 2)).Sum();

        return clonalDiversity;
    }

    public static float ComputeMeanDriversPerCell(List<SubClone> SubClones)
        => (float)SubClones.Select(clone => clone.TotalCount * clone.NumberDrivers).Sum() / SubClones.Select(clone => clone.TotalCount).Sum();

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