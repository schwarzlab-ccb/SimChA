// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.DataTypes;

namespace SimChA.Computation;

public static class TreeBuilder
{
    public static Dictionary<int, int> CreateParentMap(IEnumerable<SubClone> subClones)
        => subClones.ToDictionary(sc => sc.CloneId, sc => sc.ParentId);

    private static TreeEdge FindEdgeToParent(Dictionary<int, int> parentMap, List<SubClone> selection, int id)
    {
        int dist = 0;
        int source = id;

        do
        {
            dist++;
            source = parentMap[source];
        } while (selection.All(sc => sc.CloneId != source) && source != -1);

        return new TreeEdge { Distance = dist, SourceId = source, TargetId = id };
    }
    private static TreeEdge FindEdgeToParentWithAncestors(Dictionary<int, int> parentMap, List<SubClone> selection, List<int> internalNodes, int id)
    // Same as FindEdgeToParent but also checks whether the source is in internalNodes
    {
        int dist = 0;
        int source = id;
        do
        {
            dist++;
            source = parentMap[source];
        } while (selection.All(sc => sc.CloneId != source) && source != -1 && internalNodes.All(n => n != source));

        return new TreeEdge { Distance = dist, SourceId = source, TargetId = id };
    }
    private static List<int> FindInternalNodes(Dictionary<int, int> parentMap, List<SubClone> selection, List<SubClone> allSubClones)
    {
        Dictionary<int, int> internalNodes = new();

        foreach (var subClone in selection)
        {
            int curNode = parentMap[subClone.CloneId];
            while (selection.All(sc => sc.CloneId != curNode) && curNode != -1)
            {
                if (internalNodes.Keys.Contains(curNode))
                {
                    internalNodes[curNode]++;
                    break;
                }
                else
                {
                    internalNodes[curNode] = 0;
                    curNode = parentMap[curNode];
                }
            }
        }
        return internalNodes.Where(n => n.Value > 0 || n.Key == 0).Select(n => n.Key).ToList();

    }

    public static ParentTree BuildTree(List<SubClone> allSubClones, List<SubClone> selection)
    {
        var parentMap = TreeBuilder.CreateParentMap(allSubClones);
        List<TreeNode> nodes = new();
        List<TreeEdge> edges = new();
        int rootId = -1;

        foreach (var subClone in selection)
        {
            nodes.Add(new TreeNode { Id = subClone.CloneId, Size = subClone.AliveCount });
            edges.Add(FindEdgeToParent(parentMap, selection, subClone.CloneId));
        }

        if (edges.Count(e => e.SourceId == -1) > 1)
        {
            nodes.Add(new TreeNode { Id = -1, Size = -1 }); // Root in an abstract node since the root is missing
            rootId = -1;
        }
        else
        {
            var firstEdge = edges.Find(e => e.SourceId == -1);
            if (firstEdge != null)
            {
                edges.Remove(firstEdge);
                rootId = firstEdge.TargetId;
            }
        }

        return new ParentTree { RootId = rootId, Nodes = nodes, Edges = edges };
    }
    public static ParentTree BuildTreeWithAncestors(List<SubClone> allSubClones, float cutOff)
    {
        var parentMap = TreeBuilder.CreateParentMap(allSubClones);
        var selection = allSubClones.Where(sc => (sc.AliveCount >= cutOff)).ToList();
        var internalNodes = FindInternalNodes(parentMap, selection, allSubClones);

        List<TreeNode> nodes = new();
        List<TreeEdge> edges = new();

        foreach (var subClone in selection)
        {
            nodes.Add(new TreeNode { Id = subClone.CloneId, Size = subClone.AliveCount });
            edges.Add(FindEdgeToParentWithAncestors(parentMap, selection, internalNodes, subClone.CloneId));
        }
        foreach (var internalNode in internalNodes)
        {
            nodes.Add(new TreeNode { Id = internalNode, Size = 0 });
            edges.Add(FindEdgeToParentWithAncestors(parentMap, selection, internalNodes, internalNode));
        }

        return new ParentTree { RootId = 0, Nodes = nodes, Edges = edges.Where(e => e.TargetId != 0).ToList() };
    }
}