using SimChA.DataTypes;

namespace SimChA.Simulation;

public class Clone
{
    public Clone(int cloneId, int parentId, int distToParent, Karyotype refKaryotype, int totalMutations)
    {
        CloneId = cloneId;
        ParentId = parentId;
        DistToParent = distToParent;
        Karyotype = new Karyotype(refKaryotype);
        TotalMutations = totalMutations;
    }
    public int CloneId { get; }
    public int ParentId { get; }
    public int DistToParent { get; }
    public Karyotype Karyotype { get; set; } // This should not have a setter!
    public int TotalMutations { get; }
    public static string Header() => "clone_id\tsex\tparent\tevents\tploidy\tmissing\tmixture\tkaryotype";
    public override string ToString() 
        => $"{CloneId}\t{HGRef.Sex(Karyotype.SexXX)}\t{ParentId}\t{DistToParent}\t" +
           $"{Karyotype.CalcPloidy()}\t{Karyotype.CalcMissing()}\t{Karyotype}";
    
    public Karyotype CopyKaryotype() => new(Karyotype);

    public static (Clone root, Dictionary<int, List<Clone>> lookUp) CreateLookUp(List<Clone> clones)
    {
        Clone? root = null;
        var childLookUp = new Dictionary<int, List<Clone>>();
        foreach (var clone in clones)
        {
            if (clone.ParentId == -1 || clone.ParentId == clone.CloneId)
            {
                root = clone;
            }
            else
            {
                if (!childLookUp.ContainsKey(clone.ParentId))
                {
                    childLookUp.Add(clone.ParentId, new List<Clone>());
                }

                childLookUp[clone.ParentId].Add(clone);
            }
        }

        if (root == null)
        {
            throw new Exception("No root clone found.");
        }
        return (root, childLookUp);
    }
}