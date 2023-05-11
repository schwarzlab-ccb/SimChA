// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.EventData;
using SimChA.Simulation;

namespace SimChA.DataTypes;

public class Sample
{
    public string SampleId { get; }
    public bool SexXX { get; }
    public List<CloneIn> Clones { get; }
    public List<CNEventPars> EventPars { get;  }
    public Dictionary<int, Karyotype> Kars { get; }
    public Dictionary<int, List<CNEventDesc>> EventDescs { get; }

    public Sample(string sampleId, bool sexXX, List<CloneIn> clones,  List<CNEventPars> eventPars)
    {
        SampleId = sampleId;
        SexXX = sexXX;
        Clones = clones;
        EventPars = eventPars;
        Kars = new Dictionary<int, Karyotype>();
        EventDescs = new Dictionary<int, List<CNEventDesc>>();
    }
    
    public static string Header() => "sample_id\tsex_xx\tclone_count\tmixture";
    private string MixtureString() => EventPars.Any() ? string.Join(";", EventPars.Select(e => e.Prob)) : "-" ;
    public string ToTSV() => $"{SampleId}\t{SexXX}\t{Clones.Count}\t{MixtureString()}";
}