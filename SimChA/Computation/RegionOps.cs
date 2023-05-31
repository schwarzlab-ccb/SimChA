using SimChA.DataTypes;
using SimChA.Computation;
using SimChA.Simulation;

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

    public static List<Region> DeleteRange(List<Region> regions, long start, long end)
    {
        var newRegions = new List<Region>();
        long seekPos = 0;
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
                var newRegion = region with 
                {
                    End = region.End + start - seekPos - region.Length,
                    GeneLists = Fitness.GetGeneList(region.Start, region.End + start - seekPos - region.Length, region.ChrNo)
                };
                AddIfNotEmpty(newRegions, newRegion);
            }
            else if (start < seekPos) // start before the region, end inside the region
            {
                var newRegion = region with 
                {
                    Start = region.Start - seekPos + end,
                    GeneLists = Fitness.GetGeneList(region.Start - seekPos + end, region.End, region.ChrNo)                    
                };
                AddIfNotEmpty(newRegions, newRegion);
            }
            else // Both coordinates inside of the region
            {
                var firstRegion = region with 
                {
                    End = region.End + start - seekPos - region.Length,
                    GeneLists = Fitness.GetGeneList(region.Start, region.End + start - seekPos - region.Length, region.ChrNo)
                };
                AddIfNotEmpty(newRegions, firstRegion);

                var secondRegion = region with 
                {
                    Start = region.Start - seekPos + end,
                    GeneLists = Fitness.GetGeneList(region.Start - seekPos + end, region.End, region.ChrNo)
                };
                AddIfNotEmpty(newRegions, secondRegion);
            }

            seekPos += region.Length;
        }

        return newRegions;
    }

    public static List<Region> CopyRange(List<Region> regions, long start, long end)
    {
        var newRegions = new List<Region>();
        long seekPos = 0;
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
                var newRegion = region with 
                { 
                    Start = region.Start + start - seekPos,
                    GeneLists = Fitness.GetGeneList(region.Start + start - seekPos, region.End, region.ChrNo)
                };
                AddIfNotEmpty(newRegions, newRegion);
            }
            else if (start < seekPos) // start before the region
            {
                var newRegion = region with 
                {
                    End = region.Start + (end - seekPos),
                    GeneLists = Fitness.GetGeneList(region.Start, region.Start + (end - seekPos), region.ChrNo)
                };
                AddIfNotEmpty(newRegions, newRegion);
            }
            else // Both coordinates inside of the region
            {
                var newRegion = region with 
                {
                    Start = region.Start + start - seekPos,
                    End = region.Start + end -  seekPos,
                    GeneLists = Fitness.GetGeneList(region.Start + start - seekPos, region.Start + end - seekPos, region.ChrNo)
                }; 
                AddIfNotEmpty(newRegions, newRegion);
            }

            seekPos += region.Length;
        }

        return newRegions;
    }

    public static (List<Region>, List<Region>) SplitRegions(List<Region> regions, long pos)
    {
        long seekPos = 0;
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
                var firstPart = region with 
                { 
                    End = region.Start + pos - seekPos,
                    GeneLists = Fitness.GetGeneList(region.Start, region.Start + pos - seekPos, region.ChrNo)
                };
                AddIfNotEmpty(beforeRegions, firstPart);
                var secondPart = region with 
                { 
                    Start = firstPart.End,
                    GeneLists = Fitness.GetGeneList(firstPart.End, region.End, region.ChrNo)
                };
                AddIfNotEmpty(afterRegions, secondPart);
            }

            seekPos += region.Length;
        }

        return (beforeRegions, afterRegions);
    }

    public static List<T> StitchRegions<T>(List<T> regions) where T : GeneRegionList
    {
        var newRegions = new List<T>();
        bool[] merged = new bool[regions.Count];
        for (int i = 0; i < regions.Count; i++)
        {
            if (merged[i])
            {
                continue;
            }
            var newRegion = regions[i];
            for (int j = i + 1; j < regions.Count; j++)
            {
                if (merged[j] 
                    || regions[j].ChrNo != newRegion.ChrNo 
                    || regions[j].Start != newRegion.End)
                {
                    continue;
                }
                newRegion = newRegion with 
                {
                    End = regions[j].End,
                    GeneLists = Fitness.GetGeneList(newRegion.Start, regions[j].End, newRegion.ChrNo)
                };
                merged[j] = true;
            }
            newRegions.Add(newRegion);
        }
        return newRegions;
    }
    
    public static List<Region> GlueNeighbours(List<Region> regions)
    {
        var newRegions = new List<Region>();
        bool[] merged = new bool[regions.Count];
        for (int i = 0; i < regions.Count; i++)
        {
            if (merged[i])
            {
                continue;
            }
            var newRegion = regions[i];
            int j = i + 1; 
            if (j < regions.Count
                && !merged[j] 
                && regions[j].ChrNo == newRegion.ChrNo 
                && regions[j].Start == newRegion.End 
                && regions[j].Forward == newRegion.Forward)
            {
                newRegion = newRegion with 
                {
                    End = regions[j].End,
                    GeneLists = Fitness.GetGeneList(newRegion.Start, regions[j].End, newRegion.ChrNo)
                };
                merged[j] = true;
            }
            newRegions.Add(newRegion);
        }
        return newRegions;
    }

    public static List<Region> InvertRegions(IEnumerable<Region> regions)
        => regions.Select(r => r with { Forward = false }).Reverse().ToList();

    public static List<Region> ConcatRegions(IEnumerable<IEnumerable<Region>> listOfRegions)
        => listOfRegions.SelectMany(x => x).ToList();

    public static List<Region> ConcatRegions(IEnumerable<Region> first, IEnumerable<Region> second)
        => first.Concat(second).ToList();
    
    public static long GetLength(IEnumerable<Region> regions)
        => regions.Sum(r => r.Length);
    
    public static List<List<Region>> Scatter(List<long> locs, List<Region> regions)
    {
        if (!locs.Any()) 
        {
            return new List<List<Region>> { regions };
        }
        
        // First region
        var newRegions = new List<List<Region>> { CopyRange(regions, 0, locs[0]) };
        // Internal regions
        for (int i = 0; i < locs.Count - 1; i++)
        {
            long start = locs[i];
            long end = locs[i + 1];
            var copy = CopyRange(regions, start, end);
            newRegions.Add(copy);
        }
        // Last region
        newRegions.Add(CopyRange(regions, locs.Last(), Contig.Length(regions)));
        return newRegions;
    }
    
    public static List<Region> Gather(List<List<Region>> newRegions, IEnumerable<int> indices) 
        => ConcatRegions(indices.Select(i => newRegions[i]));
}
