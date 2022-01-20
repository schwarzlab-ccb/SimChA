using SimChA.Simulation;

namespace SimChA.DataTypes;

public class SubClone
{
    public int FirstGen;
    public int CloneId;
    public int ParentId;
    public int AliveCount => Generations.Sum(pair => pair.Value);
    public Karyotype Karyotype;
    public double DivisionRate;
    public double MutationRate;
    public Dictionary<int, int> Generations;

    public SubClone(int cloneId, int parentId, int generation, double divisionRate, double mutationRate, Karyotype karyotype) 
    {
        CloneId = cloneId;
        ParentId = parentId;
        Karyotype = karyotype;
        DivisionRate = divisionRate;
        MutationRate = mutationRate;
        FirstGen = generation;
        Generations = new Dictionary<int, int> { { FirstGen, 1 } };
    }

    public SubClone(SubClone other)
    {
        CloneId = other.CloneId;
        ParentId = other.CloneId;
        Karyotype = other.Karyotype;
        MutationRate = other.MutationRate; 
        DivisionRate = other.DivisionRate;
        Generations = other.Generations;
        FirstGen = other.FirstGen;
    }

    public SubClone CreateChild(int newId, int generation) 
        => new(this)
        {
            FirstGen = generation,
            CloneId = newId,
            ParentId = CloneId,
            Karyotype = new Karyotype(Karyotype),
            Generations = new Dictionary<int, int> { { generation, 1 } }
        };

    public override string ToString()
    {
        return $"ID:{CloneId}, Parent:{ParentId}, Cells: {AliveCount}, Karyotype: {Karyotype}";
    }
}