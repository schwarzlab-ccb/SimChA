namespace SimChA.DataTypes;

public class SubClone
{
    public int TotalCount;
    public int AliveCount;
    public int CloneId;
    public int ParentId;
    public Karyotype Karyotype;
    public float MutationRate;
    public float DivisionRate;

    public SubClone(int cloneId, int parentId, Karyotype karyotype) 
    {
        CloneId = cloneId;
        ParentId = parentId;
        Karyotype = karyotype;
        TotalCount = AliveCount = 1;
        MutationRate = DivisionRate = 1f;
    }

    public SubClone(SubClone other, int cloneId)
    {
        CloneId = cloneId;
        ParentId = other.CloneId;
        Karyotype = new Karyotype(other.Karyotype);
        TotalCount = AliveCount = 1;
        MutationRate = other.MutationRate; 
        DivisionRate = other.DivisionRate; 
    }
}