// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public class Chromosome
{
    public List<Region> Regions { get; }

    public int Length => Regions.Sum(r => r.Length);
        
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
}