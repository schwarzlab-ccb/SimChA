using SimChA.Simulation;

namespace SimChA.DataTypes;

public class SubClone
{
    public int FirstGen { get; }
    public int CloneId { get; }
    public int ParentId { get; }
    public double DivisionRate { get; }
    public int NumberDrivers { get; }
    private List<(int Alive, int Dead)> Cells { get; } 
    
    // public Karyotype Karyotype { get; }

    public int AliveCount => Cells.Last().Alive;
    public int DeadCount => Cells.Last().Dead;
    public int TotalCount => AliveCount + DeadCount;
    public int LastGen => FirstGen + Cells.Count;
    
    public SubClone(int cloneId, int parentId, int generation, double divisionRate, int numberDrivers, int popSize = 1) 
    {
        CloneId = cloneId;
        ParentId = parentId;
        NumberDrivers = numberDrivers;
        // Karyotype = new Karyotype(karyotype);
        DivisionRate = Math.Clamp(divisionRate, 0, 1);
        FirstGen = generation;
        Cells = new List<(int, int)> { (popSize, 0) };
    }
    
    public SubClone CreateChild(int newId, int generation, double divRateChange, int numberDrivers)
        => new(newId, CloneId, generation, divRateChange, numberDrivers);
    
    public override string ToString() 
        => $"ID:{CloneId}, Parent:{ParentId}, Alive: {AliveCount}, " +
           $"Dead: {DeadCount}, Drivers: {NumberDrivers}, DivisionRate: {DivisionRate}";

    public int AliveAtGen(int gen)
        => gen >= FirstGen && gen < LastGen ? Cells[gen - FirstGen].Alive : 0;

    public int TotalAtGen(int gen)
    {
        return gen >= FirstGen && gen < LastGen ? Cells[gen - FirstGen].Alive + Cells[gen - FirstGen].Dead : 0;
    }
    
    public void NewGen(int genAlive, int genDead) 
        => Cells.Add((genAlive, genDead));

    public void AddNewCells(int newAliveCount) 
        => Cells[^1] = (Cells[^1].Dead + newAliveCount, Cells[^1].Dead);
}