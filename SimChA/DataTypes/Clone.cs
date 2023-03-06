using SimChA.Simulation;

namespace SimChA.DataTypes;

public class Clone
{
    public Clone(int cloneId, int parentId, string name, int distToParent, Karyotype refKaryotype, int totalMutations)
    {
        CloneId = cloneId;
        ParentId = parentId;
        Name = name;
        DistToParent = distToParent;
        Karyotype = new Karyotype(refKaryotype);
        ChildrenIDs = new List<int>();
        TotalMutations = totalMutations;
    }
    
    public int CloneId { get; }
    public int ParentId { get; }
    public int DistToParent { get; }
    public List<int> ChildrenIDs { get; }
    public string Name { get; }
    public Karyotype Karyotype { set; get; }
    public int TotalMutations { get; }

    public override string ToString()
        => $"ID: {CloneId}, Name: {Name}, Sex: {Karyotype.Sex}, Parent: {ParentId}, Events: {DistToParent}, Coverage: {Karyotype.CalcCoverage()}, Karyotype: {Karyotype}";

    public Karyotype CopyKaryotype() => new(Karyotype);
}