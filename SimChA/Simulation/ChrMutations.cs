using SimChA.DataTypes;

namespace SimChA.Simulation;

public static class ChrMutations
{
    private static void AddIfNotEmpty(ICollection<Region> regions, Region region)
    {
        if (region.Length > 0)
        {
            regions.Add(region);
        }
    }
    
    public static List<Region> DeleteRange(List<Region> regions, int start, int end)
    {
        var newRegions = new List<Region>();
        int seekPos = 0;
        foreach (var region in regions)
        {
            if (start >= seekPos + region.Length) // region before start
            {
                AddIfNotEmpty(newRegions, region);
            }
            else if (end <= seekPos) // region after end
            {
                AddIfNotEmpty(newRegions, region);
            }
            else if (end > seekPos + region.Length) // star inside, end outside of region
            {
                var newRegion = region;
                newRegion.End = region.End + start - seekPos - region.Length;
                AddIfNotEmpty(newRegions, newRegion);
            }
            else if (start < seekPos) // start before the region, end inside the region
            {
                var newRegion = region;
                newRegion.Start = region.Start - seekPos + end;
                AddIfNotEmpty(newRegions, newRegion);
            }
            else // Both coordinates inside of the region
            {
                var firstRegion = region;
                firstRegion.End = region.End + start - seekPos - region.Length; 
                AddIfNotEmpty(newRegions, firstRegion);
                
                var secondRegion = region;
                secondRegion.Start = region.Start - seekPos + end;
                AddIfNotEmpty(newRegions, secondRegion);
            }
            seekPos += region.Length;
        }
        return newRegions;
    }

    public static List<Region> CopyRange(List<Region> regions, int start, int end)
    {
        var newRegions = new List<Region>();
        int seekPos = 0;
        foreach (var region in regions)
        {
            if (start >= seekPos + region.Length) // region before start
            {
            }
            else if (end <= seekPos) // region after end
            {
                break;
            }
            else if (end > seekPos + region.Length) // end outside of region
            {
                var newRegion = region;
                newRegion.Start = region.Start + start - seekPos;
                AddIfNotEmpty(newRegions, newRegion);
            }
            else if (start < seekPos) // start before the region
            {
                var newRegion = region;
                newRegion.End = region.End + seekPos - end + 1;
                AddIfNotEmpty(newRegions, newRegion);
            }
            else // Both coordinates inside of the region
            {
                AddIfNotEmpty(newRegions, region);
            }
            seekPos += region.Length;
        }
        return newRegions;
    }

    public static (List<Region>, List<Region>) SplitRegions(List<Region> regions, int pos)
    {
        var newRegions = new List<Region>();
        int seekPos = 0;
        var beforeRegions = new List<Region>();
        var afterRegions = new List<Region>();
        foreach (var region in regions)
        {
            if (pos > seekPos + region.Length) // region before pos
            {
                AddIfNotEmpty(beforeRegions, region);
            }
            else if (pos <= seekPos) // region after pos
            {
                AddIfNotEmpty(afterRegions, region);
            }
            else // split inside the region
            {
                var firstPart = region;
                firstPart.End -= pos - seekPos;
                AddIfNotEmpty(beforeRegions, firstPart);
                var secondPart = region;
                secondPart.Start += pos - seekPos;
                AddIfNotEmpty(afterRegions, secondPart);
            }
            seekPos += region.Length;
        }
        return (newRegions, newRegions);
    }
    
    public static List<Region> InvertRegions(List<Region> regions)
    {
        var newRegions = new List<Region>(regions);
        newRegions.Reverse();
        newRegions.ForEach(r => r.Forward = false);
        return newRegions;
    }
    
    public static List<Region> ConcatRegions(IEnumerable<List<Region>> listOfRegions)
    {
        return listOfRegions.SelectMany(x => x).ToList();
    }
}