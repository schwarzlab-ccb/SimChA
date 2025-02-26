using SimChA.Data;

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

    private static Region ResizeFront(Region region, long howMuch)
    {
        var newRegion = region with { Start = region.Start + howMuch };
        return UpdateSNVs(newRegion);
    }
    
    private static Region ResizeBack(Region region, long howMuch)
    {
        var newRegion = region with { End = region.Start + howMuch };
        return UpdateSNVs(newRegion);
    }    
    
    private static Region ResizeBoth(Region region, long start, long end)
    {
        var newRegion = region with { Start = region.Start + start, End = region.Start + end };
        return UpdateSNVs(newRegion);
    }
    
    private static Region UpdateSNVs(Region region)
    {
        if (region.SNVs == null)
        {
            return region;
        }
        var newSNVs = new List<SNV>(region.SNVs);
        var lostSNVs = region.SNVs.Where(snv => snv.Pos <= region.Start || region.End <= snv.Pos);
        foreach (var snv in lostSNVs)
        {
            newSNVs.Remove(snv);
        }
        if (newSNVs.Count == 0)
        {
            newSNVs = null;
        }
        return region with {SNVs = newSNVs};
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
                var newRegion = ResizeBack(region, start - seekPos);
                AddIfNotEmpty(newRegions, newRegion);
            }
            else if (start < seekPos) // start before the region, end inside the region
            {
                var newRegion = ResizeFront(region, end - seekPos);
                AddIfNotEmpty(newRegions, newRegion);
            }
            else // None coordinates inside of the region
            {
                var firstRegion = ResizeBack(region, start - seekPos);
                AddIfNotEmpty(newRegions, firstRegion);

                var secondRegion = ResizeFront(region, end - seekPos);
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
            return new List<Region>
            {
                Capacity = 0
            };
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
            return new List<Region>
            {
                Capacity = 0
            };
        }
        if (pArm)
        {   
            return includeCentromere ? regions.GetRange(0, index + 1) : regions.GetRange(0, index);
        }
        return includeCentromere ? regions.GetRange(index, regions.Count - index) : regions.GetRange(index + 1, regions.Count - index - 1);
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
                var newRegion = ResizeFront(region, start - seekPos);
                AddIfNotEmpty(newRegions, newRegion);
            }
            else if (start < seekPos) // start before the region
            {
                var newRegion = ResizeBack(region, end - seekPos);
                AddIfNotEmpty(newRegions, newRegion);
            }
            else // None coordinates inside of the region
            {
                var newRegion = ResizeBoth(region, start - seekPos, end - seekPos);
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
                var firstPart = ResizeBack(region, pos - seekPos);
                AddIfNotEmpty(beforeRegions, firstPart);
                var secondPart = ResizeFront(region, pos - seekPos);
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
                var newSNVs = region.SNVs ?? new List<SNV>();
                // See if the SNV is present in the region
                var snv = newSNVs.FirstOrDefault(s => s.Pos == region.Start + location - seekPos);
                if (snv != null)
                {
                    // Update the existing SNV
                    newSNVs.Remove(snv);
                    newSNVs.Add(snv with { Alt = newNucleotide });
                }
                else
                {
                    // Add a new SNV
                    newSNVs.Add(new SNV(region.Start + location - seekPos, region.Chrom, newNucleotide));
                }
                var newRegion = region with { SNVs = newSNVs };
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

    public static List<Region> InvertRegions(IEnumerable<Region> regions)
        => regions.Select(r => r with {Start = r.End * -1, End = r.Start * -1}).Reverse().ToList();

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
    
    private static List<SNV>? MergeSNVs(List<SNV>? snvs1, List<SNV>? snvs2)
    {
        if (snvs1 == null && snvs2 == null)
        {
            return null;
        }
        if (snvs1 == null)
        {
            return snvs2;
        }
        if (snvs2 == null)
        {
            return snvs1;
        }
        return snvs1.Concat(snvs2).ToList();
    }
    
    public static List<Region> MergeRegions(List<Region> regions)
    {
        var newRegions = new List<Region>();
        for (int i = 0; i < regions.Count; i++)
        {
            var cur = regions[i];
            if (i == 0)
            {
                newRegions.Add(cur);
            }
            else
            {
                var last = newRegions[^1];
                if (last.Chrom == cur.Chrom && last.End == cur.Start)
                {
                    newRegions[^1] = last with { End = cur.End, SNVs = MergeSNVs(last.SNVs, cur.SNVs) };;
                }
                else
                {
                    newRegions.Add(cur);
                }
            }
        }
        return newRegions;
    }
}
