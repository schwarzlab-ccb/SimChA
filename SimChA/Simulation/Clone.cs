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
    public static string Header() => "ID\tName\tSex\tParent\tEvents\tPloidy\tMissing\tKaryotype";
    public override string ToString() 
        => $"{CloneId}\t{Name}\t{Karyotype.Sex}\t{ParentId}\t{DistToParent}\t{Karyotype.CalcPloidy()}\t{Karyotype.CalcMissing()}\t{Karyotype}";
    public Karyotype CopyKaryotype() => new(Karyotype);
}