using SimChA.DataTypes;

namespace SimChA.Simulation;

public static class ChrMutations
{
    public static void DeleteRegion(Chromosome chr, int start, int end)
    {
        int seekPos = 0;
        for (int i = 0; i < chr.Regions.Count; i++)
        {
            var region = chr.Regions[i];
            if (start > seekPos + region.Length) // region before start
            {
                continue;
            }
            if (end > seekPos + region.Length) // end outside of region
            {
                var newRegion = region;
                newRegion.End = region.End + start - seekPos - region.Length; 
                chr.Regions[i] = newRegion;
            }
            else if (start < seekPos) // start before the region
            {
                var newRegion = region;
                newRegion.Start = region.Start - seekPos + end;
                chr.Regions[i] = newRegion;
                break;
            }
            else
            {
                var firstRegion = region;
                firstRegion.End = region.End + start - seekPos - region.Length; 
                var secondRegion = region;
                secondRegion.Start = region.Start - seekPos + end;
                chr.Regions[i] = firstRegion;
                chr.Regions.Insert(i + 1, secondRegion);
                break;
            }
            seekPos += region.Length;
        }
    }
}