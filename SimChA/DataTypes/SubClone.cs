namespace SimChA.DataTypes;

public class SubClone
{
    public SubClone(int cloneId, int parentId, int mutCount, int popSize)
    {
        CloneId = cloneId;
        ParentId = parentId;
        MutCount = mutCount;
        CellCount = popSize;
    }
    
    public int CloneId { get; }
    public int ParentId { get; }
    public int MutCount { get; }
    public int CellCount { set; get; }
    public bool IsAlive => CellCount > 0;

    public SubClone CreateChild(int newId)
        => new(newId, CloneId, MutCount + 1, 1);

    public override string ToString()
        => $"ID:{CloneId}, Parent:{ParentId}, Cells: {CellCount}, Muts: {MutCount}";
    }