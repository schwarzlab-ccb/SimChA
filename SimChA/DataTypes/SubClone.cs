using SimChA.Simulation;

namespace SimChA.DataTypes;

public class SubClone
{
    public int CloneId;
    public int ParentId;
    public int TotalCount;
    public int AliveCount;
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

    public SubClone(SubClone other)
    {
        CloneId = other.CloneId;
        ParentId = other.CloneId;
        TotalCount = other.TotalCount;
        AliveCount = other.AliveCount;
        Karyotype = other.Karyotype;
        MutationRate = other.MutationRate; 
        DivisionRate = other.DivisionRate; 
    }

    public SubClone CreateChild(int newId) 
        => new(this)
        {
            CloneId = newId,
            ParentId = CloneId,
            TotalCount = 1,
            AliveCount = 1,
            Karyotype = new Karyotype(Karyotype)
        };

    public override string ToString()
    {
        return $"ID:{CloneId}, Parent:{ParentId}, Cells: {AliveCount}, Karyotype: {Karyotype}";
    }
}