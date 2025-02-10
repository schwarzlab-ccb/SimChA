using System.Text;
using SimChA.EventData;

namespace SimChA.Data;

public class Sample
{
    public string SampleId { get; }
    public SexType Sex { get; }
    public List<CloneData> Clones { get; }
    public List<CNEventPars> EventPars { get; }
    public Dictionary<string, double> Mixture { get; }
    
    public Sample(string sampleId, SexType sex, List<CloneData> clones, List<CNEventPars>? eventPars = null, Dictionary<string, double>? mixture = null)
    {
        SampleId = sampleId;
        Sex = sex;
        Clones = clones;
        EventPars = eventPars ?? new List<CNEventPars>();
        Mixture = mixture ?? new Dictionary<string, double>();
    }
    
    public static string Header() 
        => "sample_id\tsex\tclone_count\tmixture";
    
    private string MixtureString() 
        => EventPars.Any() ? string.Join(";", Mixture.Select(pair => $"{pair.Key}:{pair.Value:f4}")) : "-";
    
    public string ToTSV() => $"{SampleId}\t" +
                             $"{Sex}\t" +
                             $"{Clones.Count}\t" +
                             $"{MixtureString()}";
    public static string HeaderAsTree() 
        => "ID\tParentID\tDistance";
    
    public string ToTSVAsTree()
    {
        var sb = new StringBuilder();
        foreach (var clone in Clones)
        {
            sb.AppendLine($"{clone.CloneId}\t{clone.ParentId}\t{clone.Distance}");
        }
        return sb.ToString();
    }
}