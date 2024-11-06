using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.Computation;

// TODO: Inverted regions currently not 
public static class RegionOps
{
    private static void AddIfNotEmpty(ICollection<Region> regions, Region region)
    {
        if (region.Length > 0)
        {
            regions.Add(region);
        }
    }

    private static Region OffsetStart(Region region, long howMuch)
    {
        var newRegion = region.Forward
            ? region with { Start = region.Start + howMuch }
            : region with { End = region.End - howMuch };
        return UpdateSNVDict(newRegion);
    }
    
    private static Region OffsetEnd(Region region, long howMuch)
    {
        var newRegion = region.Forward 
            ? region with {End = region.Start + howMuch}
            : region with {Start = region.End - howMuch};
        return UpdateSNVDict(newRegion);
    }    
    private static Region OffsetBoth(Region region, long start, long end)
    {
        var newRegion = region.Forward 
            ? region with {Start = region.Start + start, End = region.Start + end}
            : region with {Start = region.End - end, End = region.End - start };
        
        return UpdateSNVDict(newRegion);
    }
    
    private static Region UpdateSNVDict(Region region)
    {
        if (region.SNVDict == null)
        {
            return region;
        }
        var newSNVDict = new Dictionary<long, Nucleotide>(region.SNVDict);
        foreach (var snv in region.SNVDict)
        {
            if (snv.Key <= region.Start || region.End <= snv.Key)
            {
                newSNVDict.Remove(snv.Key);
            }
        }
        if (newSNVDict.Keys.Count == 0)
        {
            newSNVDict = null;
        }
        return region with {SNVDict = newSNVDict};
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
                var newRegion = OffsetEnd(region, start - seekPos);
                AddIfNotEmpty(newRegions, newRegion);
            }
            else if (start < seekPos) // start before the region, end inside the region
            {
                var newRegion = OffsetStart(region, end - seekPos);
                AddIfNotEmpty(newRegions, newRegion);
            }
            else // None coordinates inside of the region
            {
                var firstRegion = OffsetEnd(region, start - seekPos);
                AddIfNotEmpty(newRegions, firstRegion);

                var secondRegion = OffsetStart(region, end - seekPos);
                AddIfNotEmpty(newRegions, secondRegion);
            }

            seekPos += region.Length;
        }

        return newRegions;
    }

    public static List<Region> DeleteArm(List<Region> regions, int index, bool pArm, bool includeCentromere)
    {
        if (index < 0 || index >= regions.Count)
        {
            return new List<Region> { };
        }
        if (pArm)
        {
            if (includeCentromere)
            {
                regions.RemoveRange(0, index + 1);
            }
            else 
            {
                regions.RemoveRange(0, index);
            }
        }
        else
        {
            if (includeCentromere)
            {
                regions.RemoveRange(index, regions.Count - index);
            }
            else
            {
                regions.RemoveRange(index + 1, regions.Count - index - 1);
            }
        }
        return regions;
    }
    
    public static List<Region> GetArm(List<Region> regions, int index, bool pArm, bool includeCentromere)
    {
        if (index < 0 || index >= regions.Count)
        {
            return new List<Region> { };
        }
        if (pArm)
        {   
            return includeCentromere ? regions.GetRange(0, index + 1) : regions.GetRange(0, index);
        }
        else
        {
            return includeCentromere ? regions.GetRange(index, regions.Count - index) : regions.GetRange(index + 1, regions.Count - index - 1);
        }
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
                var newRegion = OffsetStart(region, start - seekPos);
                AddIfNotEmpty(newRegions, newRegion);
            }
            else if (start < seekPos) // start before the region
            {
                var newRegion = OffsetEnd(region, end - seekPos);
                AddIfNotEmpty(newRegions, newRegion);
            }
            else // None coordinates inside of the region
            {
                var newRegion = OffsetBoth(region, start - seekPos, end - seekPos);
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
                var firstPart = OffsetEnd(region, pos - seekPos);
                AddIfNotEmpty(beforeRegions, firstPart);
                var secondPart = OffsetStart(region, pos - seekPos);
                AddIfNotEmpty(afterRegions, secondPart);
            }

            seekPos += region.Length;
        }

        return (beforeRegions, afterRegions);
    }
    public static List<Region> PointMutateRegion(List<Region> regions, long location, Nucleotide newNucleotide)
    {
        long seekPos = 0;
        var newRegions = new List<Region>();
        foreach (var region in regions)
        {
            if (location >= seekPos && location < seekPos + region.Length)
            {
                var newSNVDict = region.SNVDict ?? new Dictionary<long, Nucleotide>();
                newSNVDict[region.Start + location - seekPos] = newNucleotide;
                var newRegion = region with { SNVDict = newSNVDict };
                AddIfNotEmpty(newRegions, newRegion);
            }
            else
            {
                AddIfNotEmpty(newRegions, region);
            }
            seekPos += region.Length;
        }        
        return newRegions;
    }

    public static List<Region> GlueNeighbours(List<Region> regions)
    {
        // Step 1: Sort the regions
        var sortedRegions = regions.OrderBy(r => r.ChrNo).ThenBy(r => r.Start).ThenBy(r => r.End).ToList();

        var mergedRegions = new List<Region>();
        foreach (var currentRegion in sortedRegions)
        {
            bool isMerged = false;
            for (int i = 0; i < mergedRegions.Count; i++)
            {
                var existingRegion = mergedRegions[i];
                // Check if the current region can be merged with the existing one
                if (existingRegion.ChrNo == currentRegion.ChrNo &&
                    existingRegion.Forward == currentRegion.Forward &&
                    (existingRegion.End == currentRegion.Start))
                {
                    // Merge regions by updating the existing region to encompass both
                    mergedRegions[i] = existingRegion with { End = currentRegion.End };
                    isMerged = true;
                    break;
                }
            }
            // If the current region wasn't merged, add it as a new entry
            if (!isMerged)
            {
                mergedRegions.Add(currentRegion);
            }
        }

        // Return the merged and glued regions
        return mergedRegions;
    }

    public static List<Region> InvertRegions(IEnumerable<Region> regions)
        => regions.Select(r => r with { Forward = !r.Forward }).Reverse().ToList();

    public static List<Region> ConcatRegions(IEnumerable<IEnumerable<Region>> listOfRegions)
        => listOfRegions.SelectMany(x => x).ToList();

    public static List<Region> ConcatRegions(IEnumerable<Region> first, IEnumerable<Region> second)
        => first.Concat(second).ToList();
    
    public static long GetLength(IEnumerable<Region> regions)
        => regions.Sum(r => r.Length);
    
    public static List<List<Region>> Scatter(List<long> locs, List<Region> regions)
    {
        if (locs.Count == 0) 
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

    public static (Region region, long internalLocation) FindRegion(
        List<Region> regions, long location)
    {
        long seekPos = 0;
        var region = regions[0];
        for (int i = 0; i < regions.Count; i++)
        {
            region = regions[i];
            if (location >= seekPos && location < seekPos + region.Length)
            {
                return (region, region.Start + location - seekPos);
            }
            seekPos += region.Length;
        }
        throw new Exception("Couldn't find the corresponding region of the chromsome to perform an SNV. This should not occur");
    }

    public static List<Region> MergeRegions(List<Region> regions)
    {
        var newRegions = new List<Region>();
        for (int i = 0; i < regions.Count; i++)
        {
            if (i == 0)
            {
                newRegions.Add(regions[i]);
            }
            else
            {
                var last = newRegions[^1];
                if (last.ChrNo == regions[i].ChrNo && last.End == regions[i].Start && last.Forward == regions[i].Forward)
                {
                    newRegions[^1] = last with {End = regions[i].End};
                }
                else
                {
                    newRegions.Add(regions[i]);
                }
            }
        }
        return newRegions;
    }
}
