using SimChA.Data;

namespace SimChA.Computation;

public static class CloneComp
{
    public static (CloneIn root, Dictionary<string, List<CloneIn>> lookUp) CreateLookUp(List<CloneIn> clones)
    {
        CloneIn? root = null;
        var childLookUp = new Dictionary<string, List<CloneIn>>();
        foreach (var clone in clones)
        {
            var parentId = clone.ParentId;
            if (parentId == "-1" || parentId == clone.CloneId)
            {
                root = clone;
                parentId = clone.CloneId;
            }
            if (!childLookUp.ContainsKey(parentId))
            {
                childLookUp.Add(parentId, new List<CloneIn>());
            }
            if (!childLookUp.ContainsKey(clone.CloneId))
            {
                childLookUp.Add(clone.CloneId, new List<CloneIn>());
            }
            childLookUp[parentId].Add(clone);
        }
        if (root == null)
        {
            throw new Exception("No root clone found.");
        }
        return (root, childLookUp);
    }
}