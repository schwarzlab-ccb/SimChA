// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using System.Text;
using SimChA.EventData;
using SimChA.Simulation;

namespace SimChA.DataTypes;

public class Sample
{
    public string SampleId { get; }
    public SexEnum Sex { get; }
    public List<CloneIn> Clones { get; }
    public List<CNEventPars> EventPars { get; }
    public Dictionary<string, Karyotype> Kars { get; }
    public Dictionary<string, List<CNEventDesc>> EventDescs { get; }
    public Dictionary<string, CloneStat> CloneStats { get; }
    public Dictionary<string, double> Mixture { get; }
    public Dictionary<string, Signature> Signatures { get; }
    
    public Sample(string sampleId, SexEnum sex, List<CloneIn> clones, List<CNEventPars> eventPars, 
        Dictionary<string, double> mixture, Dictionary<string, Signature> signatures)
    {
        SampleId = sampleId;
        Sex = sex;
        Clones = clones;
        EventPars = eventPars;
        Mixture = mixture;
        Signatures = signatures;
        
        Kars = new Dictionary<string, Karyotype>();
        EventDescs = new Dictionary<string, List<CNEventDesc>>();
        CloneStats = new Dictionary<string, CloneStat>();
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