using SimChA.Simulation;

namespace SimChA.DataTypes;

public class SubClone
{
    public int FirstGen;
    public int CloneId;
    public int ParentId;
    public int AliveCount => Generations.Last();
    public Karyotype Karyotype;
    public double DivisionRate;
    public double MutationRate;
    public double DriverToPassengerRate;
    public List<int> Generations;

    public SubClone(int cloneId, int parentId, int generation, double divisionRate, double mutationRate, double driverToPassengerRate, Karyotype karyotype, int popSize = 1) 
    {
        CloneId = cloneId;
        ParentId = parentId;
        Karyotype = karyotype;
        DivisionRate = divisionRate;
        MutationRate = mutationRate;
        DriverToPassengerRate = driverToPassengerRate;
        FirstGen = generation;
        Generations = new List<int> { popSize };
    }

    public SubClone(SubClone other)
    {
        CloneId = other.CloneId;
        ParentId = other.CloneId;
        Karyotype = other.Karyotype;
        MutationRate = other.MutationRate;
        DriverToPassengerRate = other.DriverToPassengerRate;
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
            Generations = new List<int> { 1 }
        };

    public override string ToString()
    {
        return $"ID:{CloneId}, Parent:{ParentId}, Cells: {AliveCount}, Karyotype: {Karyotype}";
    }
}