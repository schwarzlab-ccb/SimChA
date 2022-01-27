using SimChA.Simulation;

namespace SimChA.DataTypes;

public class SubClone
{
    public int FirstGen { get; private init; }
    public int CloneId { get; private init; }
    public int ParentId { get; private init; }
    public Karyotype Karyotype { get; private init; }
    public double DivisionRate { get; private init; }
    public List<int> AliveCells { get; private init; }
    public List<int> DeadCells { get; private init; }
    
    
    public int AliveCount => AliveCells.Last();
    public int DeadCount => DeadCells.Last();
    public int TotalCount => AliveCount + DeadCount;

    public SubClone(int cloneId, int parentId, int generation, double divisionRate, Karyotype karyotype, int popSize = 1) 
    {
        CloneId = cloneId;
        ParentId = parentId;
        Karyotype = karyotype;
        DivisionRate = divisionRate;
        FirstGen = generation;
        AliveCells = new List<int> { popSize };
        DeadCells = new List<int> { 0 };
    }

    public SubClone(SubClone other)
    {
        CloneId = other.CloneId;
        ParentId = other.CloneId;
        Karyotype = other.Karyotype;
        DivisionRate = other.DivisionRate;
        AliveCells = other.AliveCells;
        DeadCells = other.DeadCells;
        FirstGen = other.FirstGen;
    }

    public SubClone CreateChild(int newId, int generation, double divRateChange) 
        => new(this)
        {
            FirstGen = generation,
            CloneId = newId,
            ParentId = CloneId,
            Karyotype = new Karyotype(Karyotype),
            AliveCells = new List<int> { 1 },
            DeadCells = new List<int> { 0 },
            DivisionRate = Math.Clamp(DivisionRate * divRateChange, 0, 1)
        };

    public override string ToString()
    {
        return $"ID:{CloneId}, Parent:{ParentId}, Alive: {AliveCount}, Dead: {DeadCount}, Karyotype: {Karyotype}";
    }

    public int MaxPopulation()
        => AliveCells.Max();

    public int PopAtGeneration(int gen) 
        => gen < FirstGen || gen >= FirstGen + AliveCells.Count ? -1 : AliveCells[gen - FirstGen];
}