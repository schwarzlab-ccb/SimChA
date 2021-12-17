namespace SimChA.DataTypes;

public class SubClone
{
    public int TotalCount;
    public int AliveCount;
    public int CloneId;
    public int ParentId;
    public Karyotype Karyotype;
    public double DivisionRate;
    public double MutationRate;

    public SubClone(int cloneId, int parentId,  double divisionRate,double mutationRate, Karyotype karyotype) 
    {
        CloneId = cloneId;
        ParentId = parentId;
        Karyotype = karyotype;
        DivisionRate = divisionRate;
        MutationRate = mutationRate;
        TotalCount = AliveCount = 1;
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