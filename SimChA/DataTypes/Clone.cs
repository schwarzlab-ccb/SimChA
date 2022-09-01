using SimChA.Simulation;

namespace SimChA.DataTypes;

public class Clone
{
    public Clone(int cloneId, int parentId, int mutCount, int popSize, Karyotype? refKaryotype)
    {
        CloneId = cloneId;
        ParentId = parentId;
        MutCount = mutCount;
        CellCount = popSize;
        Karyotype = new Karyotype(refKaryotype);
    }

    public int CloneId { get; }
    public int ParentId { get; }
    public int MutCount { get; }
    public int CellCount { set; get; }
    public bool IsAlive => CellCount > 0;
    public Karyotype? Karyotype { set; get; }

    public Clone CreateChild(int newId)
        => new(newId, CloneId, MutCount + 1, 1, Karyotype);

    public override string ToString()
        => $"ID:{CloneId}, Parent:{ParentId}, Cells: {CellCount}, Muts: {MutCount}, Karyotype: {Karyotype}";


    public Karyotype SetKaryotype() => new(Karyotype);
}