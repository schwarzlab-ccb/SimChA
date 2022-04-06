using SimChA.Simulation;

namespace SimChA.DataTypes;

public class SubClone
{
    public int FirstGen { get; }
    public int CloneId { get; }
    public int ParentId { get; }
    public double DivisionRate { get; }
    public int NumberDrivers { get; }
    private List<(long Alive, long Dead, long Decayed)> Cells { get; }
    
    public long AliveCount => Cells.Last().Alive;
    public long DeadCount => Cells.Last().Dead;
    public long DecayedCount => Cells.Last().Decayed;
    public long SampleCount => AliveCount + DeadCount;
    public long TotalCount => AliveCount + DeadCount + DecayedCount;
    public int LastGen => FirstGen + Cells.Count;
    
    public SubClone(int cloneId, int parentId, int generation, double divisionRate, int numberDrivers = 1, uint popSize = 1) 
    {
        CloneId = cloneId;
        ParentId = parentId;
        NumberDrivers = numberDrivers;
        DivisionRate = Math.Clamp(divisionRate, 0, 1);
        FirstGen = generation;
        Cells = new List<(long, long, long)> { (popSize, 0, 0) };
    }
    
    public SubClone CreateChild(int newId, int generation, double divRateChange, int numberDrivers)
        => new(newId, CloneId, generation, divRateChange, numberDrivers);
    
    public override string ToString() 
        => $"ID:{CloneId}, Parent:{ParentId}, Alive: {AliveCount}, " +
           $"Dead: {DeadCount}, Drivers: {NumberDrivers}, DivisionRate: {DivisionRate}";

    public long SampleAtGen(int gen)
        => gen >= FirstGen && gen < LastGen ? Cells[gen - FirstGen].Alive + Cells[gen - FirstGen].Dead : 0;
    
    public long AliveAtGen(int gen)
        => gen >= FirstGen && gen < LastGen ? Cells[gen - FirstGen].Alive : 0;

    public long TotalAtGen(int gen)
    {
        if (gen < FirstGen)
        {
            return 0;
        }
        if (gen >= LastGen)
        {
            return Cells.Last().Alive + Cells.Last().Dead;
        }

        return Cells[gen - FirstGen].Alive + Cells[gen - FirstGen].Dead;
    }
    
    public void NewGen(uint genAlive, uint genDead, uint decayed) 
        => Cells.Add((genAlive, genDead, decayed));
}