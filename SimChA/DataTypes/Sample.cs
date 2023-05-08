// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.Simulation;

namespace SimChA.DataTypes;

public record Sample(string SampleId, bool SexXX, List<Clone> Clones, List<CNEventP>? EventPs = null)
{
    public static string Header() => "sample_id\tsex_xx\tclone_count\tmixture";
    private string MixtureString() => EventPs == null ? "-" : string.Join(";", EventPs.Select(e => e.Prob));
    public string ToTsv() => $"{SampleId}\t{SexXX}\t{Clones.Count}\t{MixtureString()}";
}