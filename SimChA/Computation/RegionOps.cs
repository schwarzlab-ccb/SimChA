using SimChA.DataTypes;

namespace SimChA.Computation;

public static class RegionOps
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
            else if (start < seekPos && end > seekPos + region.Length) // Whole region inside deletion
            {
            }
            else if (end > seekPos + region.Length) // star inside, end outside of region
            {
                var newRegion = region with {End = region.End + start - seekPos - region.Length};
                AddIfNotEmpty(newRegions, newRegion);
            }
            else if (start < seekPos) // start before the region, end inside the region
            {
                var newRegion = region with {Start = region.Start - seekPos + end};
                AddIfNotEmpty(newRegions, newRegion);
            }
            else // Both coordinates inside of the region
            {
                var firstRegion = region with {End = region.End + start - seekPos - region.Length};
                AddIfNotEmpty(newRegions, firstRegion);

                var secondRegion = region with {Start = region.Start - seekPos + end};
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
            else if (start < seekPos && end > seekPos + region.Length) // Whole region inside copy
            {
                AddIfNotEmpty(newRegions, region);
            }
            else if (end > seekPos + region.Length) // end outside of region
            {
                var newRegion = region with { Start = region.Start + start - seekPos};
                AddIfNotEmpty(newRegions, newRegion);
            }
            else if (start < seekPos) // start before the region
            {
                var newRegion = region with {End = region.Start + (end - seekPos)};
                AddIfNotEmpty(newRegions, newRegion);
            }
            else // Both coordinates inside of the region
            {
                var newRegion = region with {Start = region.Start + start - seekPos, End = region.Start + end -  seekPos}; 
                AddIfNotEmpty(newRegions, newRegion);
            }

            seekPos += region.Length;
        }

        return newRegions;
    }

    public static (List<Region>, List<Region>) SplitRegions(List<Region> regions, int pos)
    {
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
                var firstPart = region with { End = region.Start + pos - seekPos };
                AddIfNotEmpty(beforeRegions, firstPart);
                var secondPart = region with { Start = firstPart.End};
                AddIfNotEmpty(afterRegions, secondPart);
            }

            seekPos += region.Length;
        }

        return (beforeRegions, afterRegions);
    }

    public static List<Region> InvertRegions(IEnumerable<Region> regions)
        => regions.Select(r => r with { Forward = false }).Reverse().ToList();

    public static List<Region> ConcatRegions(IEnumerable<IEnumerable<Region>> listOfRegions)
        => listOfRegions.SelectMany(x => x).ToList();

    public static List<Region> ConcatRegions(IEnumerable<Region> first, IEnumerable<Region> second)
        => first.Concat(second).ToList();
    
    public static long GetLength(IEnumerable<Region> regions)
        => regions.Sum(r => r.Length);
}