using SimChA.Simulation;

namespace SimChA.DataTypes;

public class SubClone
{
    public int FirstGen { get; }
    public int CloneId { get; }
    public int ParentId { get; }
    public double BirthRate { get; }
    public double DeathRate { get; }
    public int NumberDrivers { get; }
    private List<(long Alive, long Necro)> Cells { get; }
    
    public long AliveCount => Cells.Last().Alive;
    public long NecroCount => Cells.Last().Necro;
    public int LastGen => FirstGen + Cells.Count;
    public long LostCount { get; private set; }
    public long TotalCount => AliveCount + NecroCount + LostCount;
    
    public SubClone(int cloneId, int parentId, int generation, double birthRate, double deathRate, int numberDrivers = 1, int popSize = 1) 
    {
        CloneId = cloneId;
        ParentId = parentId;
        NumberDrivers = numberDrivers;
        BirthRate = birthRate;
        DeathRate = deathRate;
        FirstGen = generation;
        Cells = new List<(long, long)> { (popSize, 0) };
    }
    
    public SubClone CreateChild(int newId, int generation, double divRateChange, double deathChange, int numberDrivers)
        => new(newId, CloneId, generation, divRateChange, deathChange, numberDrivers);
    
    public override string ToString() 
        => $"ID:{CloneId}, Parent:{ParentId}, Alive: {AliveCount}, Necrotic: {NecroCount}, Lost: {LostCount}, " +
           $"Drivers: {NumberDrivers}, DivisionRate: {BirthRate}, DeathRate: {DeathRate}";
    
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
            return Cells.Last().Alive + Cells.Last().Necro;
        }

        return Cells[gen - FirstGen].Alive + Cells[gen - FirstGen].Necro;
    }

    public void NewGen(uint genAlive, uint genDead, uint genDis)
    {
        Cells.Add((genAlive, genDead));
        LostCount += genDis;
    }
}