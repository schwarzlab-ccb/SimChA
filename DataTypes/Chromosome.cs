// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public class Chromosome
{
    private readonly List<Region> _regions;

    public int Length => _regions.Sum(r => r.Length);
        
    public Chromosome(Region initialRegion)
    {
        _regions = new List<Region> { initialRegion };
    }
    
    public Chromosome(Chromosome other)
    {
        _regions = new List<Region>(other._regions);
    }

    public void DeleteRegion(int start, int end)
    {
        int seekPos = 0;
        for (int i = 0; i < _regions.Count; i++)
        {
            var region = _regions[i];
            if (start > seekPos + region.Length) // region before start
            {
                continue;
            }
            if (end > seekPos + region.Length) // end outside of region
            {
                var newRegion = region;
                newRegion.End = region.End + start - seekPos - region.Length; 
                _regions[i] = newRegion;
            }
            else if (start < seekPos) // start before the region
            {
                var newRegion = region;
                newRegion.Start = region.Start - seekPos + end;
                _regions[i] = newRegion;
                break;
            }
            else
            {
                var firstRegion = region;
                firstRegion.End = region.End + start - seekPos - region.Length; 
                var secondRegion = region;
                secondRegion.Start = region.Start - seekPos + end;
                _regions[i] = firstRegion;
                _regions.Insert(i + 1, secondRegion);
                break;
            }
            seekPos += region.Length;
        }
    }
}