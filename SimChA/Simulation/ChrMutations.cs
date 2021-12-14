using SimChA.DataTypes;

namespace SimChA.Simulation;

public static class ChrMutations
{
    public static Chromosome DeleteRegion(Chromosome chr, int start, int end)
    {
        var newRegions = new List<Region>();
        int seekPos = 0;
        foreach (var region in chr.Regions)
        {
            if (start > seekPos + region.Length) // region before start
            {
                newRegions.Add(region);
            }
            else if (end > seekPos + region.Length) // end outside of region
            {
                var newRegion = region;
                newRegion.End = region.End + start - seekPos - region.Length;
                newRegions.Add(newRegion);
            }
            else if (start < seekPos) // start before the region
            {
                var newRegion = region;
                newRegion.Start = region.Start - seekPos + end;
                newRegions.Add(newRegion);
            }
            else
            {
                var firstRegion = region;
                firstRegion.End = region.End + start - seekPos - region.Length; 
                var secondRegion = region;
                secondRegion.Start = region.Start - seekPos + end;
                newRegions.Add(firstRegion);
                newRegions.Add(secondRegion);
            }
            seekPos += region.Length;
        }
        return new Chromosome(newRegions);
    }
}