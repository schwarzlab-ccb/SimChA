using SimChA.Simulation;

namespace SimChA.DataTypes;

public class Clone
{
    public Clone(int cloneId, int parentId, string name, int distToParent, Karyotype refKaryotype)
    {
        CloneId = cloneId;
        ParentId = parentId;
        Name = name;
        DistToParent = distToParent;
        Karyotype = new Karyotype(refKaryotype);
        ChildrenIDs = new List<int>();
    }

    public int CloneId { get; }
    public int ParentId { get; }
    public int DistToParent { get; }
    public Karyotype Karyotype { set; get; }
    public List<int> ChildrenIDs { set; get; }
    public string Name { set; get; }
    public float Fitness { get; set;}

    public override string ToString()
        => $"ID:{CloneId}, Name:{Name} Parent:{ParentId},  Muts: {DistToParent}, Karyotype: {Karyotype}";

    public Karyotype CopyKaryotype() => new(Karyotype);
}