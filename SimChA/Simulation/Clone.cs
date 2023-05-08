using SimChA.DataTypes;

namespace SimChA.Simulation;

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
    public static string Header() => "clone_id\tname\tsex\tparent\tevents\tploidy\tmissing\tmixture\tkaryotype";
    public override string ToString() 
        => $"{CloneId}\t{Name}\t{HGRef.Sex(Karyotype.SexXX)}\t{ParentId}\t{DistToParent}\t" +
           $"{Karyotype.CalcPloidy()}\t{Karyotype.CalcMissing()}\t{Karyotype}";
    public Karyotype CopyKaryotype() => new(Karyotype);
}