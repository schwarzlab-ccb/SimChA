using SimChA.Computation;

namespace SimChA.Data;

public class PresentGenes
{
    public Dictionary<GeneListType, List<Gene>> Genes;

    public PresentGenes()
    {
        Genes =  GetEmptyGenes();
    }

    public PresentGenes(Dictionary<GeneListType, List<Gene>> presentGenes)
        => Genes = new Dictionary<GeneListType, List<Gene>>(presentGenes
            .ToDictionary(
                kvp => kvp.Key, 
                kvp => new List<Gene>(kvp.Value)
            ));

    public PresentGenes(PresentGenes other)
        => Genes = new Dictionary<GeneListType, List<Gene>>(other.Genes
            .ToDictionary(
                kvp => kvp.Key, 
                kvp => new List<Gene>(kvp.Value)
            ));

    public static Dictionary<GeneListType, List<Gene>> GetEmptyGenes()
    {
        return new Dictionary<GeneListType, List<Gene>> {
                { GeneListType.TumorSuppressor, new List<Gene>()},
                { GeneListType.Oncogene, new List<Gene>()},
                { GeneListType.Essentiality, new List<Gene>()},
            };
    }

    public static PresentGenes CollectGenes(List<Region> regions)
    {
        var dicts = regions.Select(r => r.PresentGenes.Genes);
        return new PresentGenes(dicts
            .SelectMany(dict => dict) // Flatten all key-value pairs
            .GroupBy(kvp => kvp.Key)  // Group by key
            .ToDictionary(
                g => g.Key,
                g => g.SelectMany(kvp => kvp.Value).ToList() // Merge lists
            ));
    }
    public static PresentGenes CollectGenes(List<Contig> contigs)
    {
        var dicts = contigs.Select(c => c.PresentGenes.Genes);
        return new PresentGenes(dicts
            .SelectMany(dict => dict) // Flatten all key-value pairs
            .GroupBy(kvp => kvp.Key)  // Group by key
            .ToDictionary(
                g => g.Key,
                g => g.SelectMany(kvp => kvp.Value).ToList() // Merge lists
            ));
    }

    public static Dictionary<GeneListType, Dictionary<Gene, int>> GetGeneCounts(List<Contig> contigs)
    {
        return Enum.GetValues<GeneListType>()
            .ToDictionary(
                type => type,
                type => contigs
                    .SelectMany(c => c.PresentGenes.Genes.GetValueOrDefault(type, []) ?? []) // Get the List<Gene>
                    .GroupBy(gene => gene) // Group by gene
                    .ToDictionary(g => g.Key, g => g.Count()) // Count occurrences
            );
    }
}