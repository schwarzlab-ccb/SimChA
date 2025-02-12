// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

namespace SimChA.Data;

public record SampleStats(
    string SampleId,
    double Ploidy,
    double Coverage,
    double Fitness,
    double FitnessTarget,
    double Stress,
    double Tsg,
    double Og,
    double Ess,
    int MutCount,
    double Hemizygosity,
    double Nullizygosity)
{
    public static string Header() 
        => "sample_id" +
           "\tploidy" +
           "\tcoverage" +
           "\tfitness" +
           "\tfitness_target" +
           "\tstress" +
           "\ttsg" +
           "\tog" +
           "\tess" +
           "\tdist" +
           "\themizygosity" +
           "\tnullizygosity";
    
    public override string ToString() =>
        $"{SampleId}" +
        $"\t{Ploidy}" +
        $"\t{Coverage}" +
        $"\t{Fitness}" +
        $"\t{FitnessTarget}" +
        $"\t{Stress}" +
        $"\t{Tsg}" +
        $"\t{Og}" +
        $"\t{Ess}" +
        $"\t{MutCount}" +
        $"\t{Hemizygosity}" +
        $"\t{Nullizygosity}";
}