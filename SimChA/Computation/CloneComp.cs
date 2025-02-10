using SimChA.Data;

namespace SimChA.Computation;

public static class CloneComp
{
    public static (CloneData root, Dictionary<string, List<CloneData>> lookUp) CreateLookUp(List<CloneData> clones)
    {
        CloneData? root = null;
        var childLookUp = new Dictionary<string, List<CloneData>>();
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
                childLookUp.Add(parentId, new List<CloneData>());
            }
            if (!childLookUp.ContainsKey(clone.CloneId))
            {
                childLookUp.Add(clone.CloneId, new List<CloneData>());
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