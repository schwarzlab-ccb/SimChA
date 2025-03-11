namespace SimChA.Data;

public class PresentGenes
{
    public Dictionary<GeneListType, List<Gene>> Genes { get; }

    public PresentGenes() =>
        Genes = new Dictionary<GeneListType, List<Gene>>
        {
            {GeneListType.TumorSuppressor, []},
            {GeneListType.Oncogene, []},
            {GeneListType.Essentiality, []},
        };

    public PresentGenes(Dictionary<GeneListType, List<Gene>> presentGenes)
        => Genes = presentGenes.ToDictionary(
            kvp => kvp.Key,
            kvp => new List<Gene>(kvp.Value)
        );

    public PresentGenes(PresentGenes other)
        => Genes = other.Genes.ToDictionary(
            kvp => kvp.Key,
            kvp => new List<Gene>(kvp.Value)
        );

    public static PresentGenes CollectGenes(List<Region> regions)
        => new(regions
            .SelectMany(r => r.PresentGenes.Genes)
            .ToLookup(kvp => kvp.Key, kvp => kvp.Value)
            .ToDictionary(g => g.Key, g => g.SelectMany(l => l).ToList()));

    public static Dictionary<GeneListType, Dictionary<Gene, int>> GetGeneCounts(List<Contig> contigs)
        => Enum.GetValues<GeneListType>()
            .ToDictionary(
                type => type,
                type => contigs
                    .SelectMany(c => c.PresentGenes.Genes.GetValueOrDefault(type, []) ?? [])
                    .ToLookup(gene => gene)
                    .ToDictionary(g => g.Key, g => g.Count()));
    
    public static void UpdateGeneCounts(
        Dictionary<GeneListType, Dictionary<Gene, int>> current,
        Dictionary<GeneListType, List<Gene>>? toRemove,
        Dictionary<GeneListType, List<Gene>>? toAdd
    )
    {
        foreach (var (type, geneDict) in current)
        {
            if (toRemove != null)
            {
                foreach (var gene in toRemove.GetValueOrDefault(type, []))
                {
                    geneDict[gene] -= 1;
                }
            }
            if (toAdd != null)
            {
                foreach (var gene in toAdd.GetValueOrDefault(type, []))
                {
                    geneDict[gene] += 1;
                }
            }
        }
    }

    public static void DoubleGeneCounts(Dictionary<GeneListType, Dictionary<Gene, int>> gc)
    {
        foreach (var (_, geneCounts) in gc)
        {
            foreach (var gene in geneCounts.Keys)
            {
                geneCounts[gene] *= 2;
            }
        }
    }
}