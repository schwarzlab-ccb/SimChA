using SimChA.DataTypes;

namespace SimChA.Simulation;

public record CloneIn(int CloneId, int ParentId, int Distance, double FitnessTarget)
{
    public static string Header() => "clone_id\tsex\tparent\tdistance\tploidy\tmissing\tkaryotype";
    public override string ToString()
        => $"{CloneId}\t{ParentId}\t{Distance}\t";
    // + $"{Karyotype.CalcPloidy()}\t{Karyotype.CalcMissing()}\t{Karyotype}";

    public static (CloneIn root, Dictionary<int, List<CloneIn>> lookUp) CreateLookUp(List<CloneIn> clones)
    {
        CloneIn? root = null;
        var childLookUp = new Dictionary<int, List<CloneIn>>();
        foreach (var clone in clones)
        {
            int parentId = clone.ParentId;
            if (parentId == -1 || parentId == clone.CloneId)
            {
                root = clone;
                parentId = clone.CloneId;
            }
            if (!childLookUp.ContainsKey(parentId))
            {
                childLookUp.Add(parentId, new List<CloneIn>());
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