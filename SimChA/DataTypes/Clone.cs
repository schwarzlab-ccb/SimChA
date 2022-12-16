using SimChA.Simulation;

namespace SimChA.DataTypes;

public class Clone
{
    public Clone(int cloneId, int parentId, string name, int mutCount, Karyotype? refKaryotype)
    {
        CloneId = cloneId;
        ParentId = parentId;
        Name = name;
        MutCount = mutCount;
        Karyotype = new Karyotype(refKaryotype);
        ChildrenIDs = new List<int>();
    }

    public int CloneId { get; }
    public int ParentId { get; }
    public int MutCount { get; }
    public Karyotype? Karyotype { set; get; }
    public List<int> ChildrenIDs { set; get; }
    public string Name { set; get; }
    public float? deltaFitness { get; set;}

    public override string ToString()
        => $"ID:{CloneId}, Name:{Name} Parent:{ParentId},  Muts: {MutCount}, Karyotype: {Karyotype}";

    public Karyotype SetKaryotype() => new(Karyotype);
}