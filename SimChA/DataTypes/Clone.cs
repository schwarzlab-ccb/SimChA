using SimChA.Simulation;

namespace SimChA.DataTypes;

public class Clone
{
    public Clone(int cloneId, int parentId, string name, int mutCount, int popSize, Karyotype? refKaryotype)
    {
        CloneId = cloneId;
        ParentId = parentId;
        Name = name;
        MutCount = mutCount;
        CellCount = popSize;
        Karyotype = new Karyotype(refKaryotype);
        ChildrenIDs = new List<int>();
    }

    public int CloneId { get; }
    public int ParentId { get; }
    public int MutCount { get; }
    public int CellCount { set; get; }
    public bool IsAlive => CellCount > 0;
    public Karyotype? Karyotype { set; get; }
    public List<int> ChildrenIDs { set; get; }
    public string Name { set; get; }

    public Clone CreateChild(int newId)
        => new(newId, CloneId, newId.ToString(), MutCount + 1, 1, Karyotype);

    public override string ToString()
        => $"ID:{CloneId}, Parent:{ParentId}, Cells: {CellCount}, Muts: {MutCount}, Karyotype: {Karyotype}";


    public Karyotype SetKaryotype() => new(Karyotype);
}