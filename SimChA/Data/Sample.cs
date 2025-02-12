using SimChA.EventData;

namespace SimChA.Data;

public class Sample
{
    public string SampleId { get; }
    public string ParentId { get; }
    public Karyotype Karyotype { get; }
    public List<CNEventDesc> Events { get; }
    public Dictionary<string, double> Mixture { get; }

    public Sample(string sampleId, string parentId, Karyotype kar, List<CNEventDesc>? events = null, Dictionary<string, double>? mixture = null)
    {
        SampleId = sampleId;
        ParentId = parentId;
        Karyotype = kar;
        Events = events ?? new List<CNEventDesc>();
        Mixture = mixture ?? new Dictionary<string, double>();
    }
    
    public static string Header() 
        => "sample_id\tparent_id\tsex\tdistance\tfitness\tmixture";
    
    public string ToTSV() 
        => $"{SampleId}\t{ParentId}\t{Karyotype.Sex}\t{Events.Count}\t{Karyotype.FitnessVal:f4}\t{MixtureString()}";
     
    private string MixtureString() 
        => string.Join(";", Mixture.Select(pair => $"{pair.Key}:{pair.Value:f4}"));
}