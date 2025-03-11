namespace SimChA.Data;

public class PresentGenes
{
    public Dictionary<GeneListType, List<Gene>> Genes { get;  }

    public PresentGenes() =>
        Genes = new Dictionary<GeneListType, List<Gene>>
        {
            { GeneListType.TumorSuppressor, []},
            { GeneListType.Oncogene, []},
            { GeneListType.Essentiality, []},
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

    public static Dictionary<GeneListType, Dictionary<Gene, int>> UpdateGeneCounts(
        Dictionary<GeneListType, Dictionary<Gene, int>> gc,
        Dictionary<GeneListType, List<Gene>> genesBefore,
        Dictionary<GeneListType, List<Gene>> genesAfter
    )
    {
        var d = gc.ToDictionary(kvp => kvp.Key,
            kvp => 
            {
                var type = kvp.Key;
                var countBefore = genesBefore.TryGetValue(type, out var genes1) 
                    ? genes1.GroupBy(g => g)
                            .ToDictionary(g => g.Key, g => g.Count()) 
                    : new Dictionary<Gene, int>();
                var countAfter = genesAfter.TryGetValue(type, out var genes2) 
                    ? genes2.GroupBy(g => g)
                            .ToDictionary(g => g.Key, g => g.Count()) 
                    : new Dictionary<Gene, int>();

                 // Collect all genes that appear in either countBefore or countAfter
                var genesToUpdate = new HashSet<Gene>(countBefore.Keys);
                genesToUpdate.UnionWith(countAfter.Keys);

                var updatedDict = kvp.Value.ToDictionary(geneCN => geneCN.Key, geneCN => geneCN.Value);

                foreach (var g in genesToUpdate)
                {
                    updatedDict[g] = kvp.Value.GetValueOrDefault(g, 0) 
                                + countAfter.GetValueOrDefault(g, 0) 
                                - countBefore.GetValueOrDefault(g, 0);
                }
                return updatedDict;
            });
        return d;
    }
    
    public static Dictionary<GeneListType, Dictionary<Gene, int>> UpdateGeneCounts(
        Dictionary<GeneListType, Dictionary<Gene, int>> gc,
        Dictionary<GeneListType, List<Gene>> genes, bool subtract)
    {
        int sign = subtract ? -1 : 1;
        var d = gc.ToDictionary(kvp => kvp.Key,
            kvp => 
            {
                var type = kvp.Key;
                var countBefore = genes.TryGetValue(type, out var gL) 
                    ? gL.GroupBy(g => g)
                            .ToDictionary(g => g.Key, g => g.Count()) 
                    : new Dictionary<Gene, int>();
                 // Collect all genes that appear in either countBefore or countAfter
                var genesToUpdate = new HashSet<Gene>(countBefore.Keys);
                
                var updatedDict = kvp.Value.ToDictionary(geneCN => geneCN.Key, geneCN => geneCN.Value);

                foreach (var g in genesToUpdate)
                {
                    updatedDict[g] = kvp.Value.GetValueOrDefault(g, 0) 
                                + sign * countBefore.GetValueOrDefault(g, 0);
                }
                return updatedDict;
            });
        return d;
    }

    public static  Dictionary<GeneListType, Dictionary<Gene, int>> DoubleGeneCounts(
        Dictionary<GeneListType, Dictionary<Gene, int>> gc)
    {
        var d = gc.ToDictionary(
            kvp => kvp.Key,
            kvp => kvp.Value
                .ToDictionary(geneCN => geneCN.Key, geneCN => geneCN.Value * 2)
            );
        return d;
    }
    
}