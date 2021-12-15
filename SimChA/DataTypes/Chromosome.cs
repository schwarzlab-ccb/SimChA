// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.Simulation;

namespace SimChA.DataTypes;

public class Chromosome
{
    private List<Region> Regions { get; set; }

    public int Length => Regions.Sum(r => r.Length);

    public Chromosome(IEnumerable<Region> regions)
    {
        Regions = regions.Where(r => r.Length > 0).ToList();
    }
    
    public Chromosome(Region initialRegion)
    {
        Regions = new List<Region> { initialRegion };
    }
    
    public Chromosome(Chromosome other)
    {
        Regions = new List<Region>(other.Regions);
    }

    public override string ToString() 
        => "[" + string.Join(",", Regions.Select(r => r.ToString())) + "]";

    public void DeleteRange(int start, int end)
    {
        Regions = ChrMutations.DeleteRange(Regions, start, end);
    }
}