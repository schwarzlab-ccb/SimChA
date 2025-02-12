using SimChA.EventData;

namespace SimChA.Data;

public class Sample
{
    public string SampleId { get; }
    public string ParentId { get; }
    public Karyotype Karyotype { get; }
    public List<CNEventDesc> Events { get; }
    private Dictionary<string, double> Mixture { get; }

    public Sample(string sampleId, string parentId, Karyotype kar, List<CNEventDesc>? events = null, Dictionary<string, double>? mixture = null)
    {
        SampleId = sampleId;
        ParentId = parentId;
        Karyotype = kar;
        Events = events ?? new List<CNEventDesc>();
        Mixture = mixture ?? new Dictionary<string, double>();
    }
    
    public string MixtureString() 
        => string.Join(";", Mixture.Select(pair => $"{pair.Key}:{pair.Value:f4}"));
}