using SimChA.Simulation;

namespace SimChA.DataTypes;

public class SubClone
{
    public int FirstGen { get; private init; }
    public int CloneId { get; private init; }
    public int ParentId { get; private init; }
    public Karyotype Karyotype { get; private init; }
    public double DivisionRate { get; private init; }
    private List<(int, int)> Cells { get; init; } // (Alive,Dead)

    public int AliveCount => Cells.Last().Item1;
    public int DeadCount => Cells.Last().Item2;
    public int TotalCount => AliveCount + DeadCount;
    public int LastGen => FirstGen + Cells.Count;
    
    public SubClone(int cloneId, int parentId, int generation, double divisionRate, Karyotype karyotype, int popSize = 1) 
    {
        CloneId = cloneId;
        ParentId = parentId;
        Karyotype = karyotype;
        DivisionRate = Math.Clamp(divisionRate, 0, 1);
        FirstGen = generation;
        Cells = new List<(int, int)> { (popSize, 0) };
    }
    
    public SubClone CreateChild(int newId, int generation, double divRateChange)
        => new(newId, CloneId, generation, DivisionRate * divRateChange, new Karyotype(Karyotype));
    
    public override string ToString() 
        => $"ID:{CloneId}, Parent:{ParentId}, Alive: {AliveCount}, Dead: {DeadCount}, Karyotype: {Karyotype}";

    public int AliveAtGen(int gen)
        => gen >= FirstGen && gen < LastGen ? Cells[gen - FirstGen].Item1 : 0;

    public void NewGen(int genAlive, int genDead) 
        => Cells.Add((genAlive, genDead));

    public void AddNewCells(int newAliveCount) 
        => Cells[^1] = (Cells[^1].Item1 + newAliveCount, Cells[^1].Item2);
}