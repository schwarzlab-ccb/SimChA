using SimChA.IO;
using SimChA.Computation;
namespace SimChA.Data;

public record SampleStat(
    string SampleId,
    string ParentId,
    SexType Sex,
    double Ploidy,
    double FitnessVal,
    double FitnessTarget,
    double Stress,
    double Tsg,
    double Og,
    double Ess,
    int MutCount,
    double Hemizygosity,
    double Nullizygosity,
    string Mixture)
{
    public static string Header() 
        => "sample_id" +
           "\tparent_id" +
           "\tsex" +
           "\tploidy" +
           "\tfitness" +
           "\tfitness_target" +
           "\tstress" +
           "\ttsg" +
           "\tog" +
           "\tess" +
           "\tdist" +
           "\themizygosity" +
           "\tnullizygosity" + 
           "\tmixture";
    
    public override string ToString() =>
        $"{SampleId}" +
        $"\t{ParentId}" +
        $"\t{Sex}" +
        $"\t{Ploidy}" +
        $"\t{FitnessVal}" +
        $"\t{FitnessTarget}" +
        $"\t{Stress}" +
        $"\t{Tsg}" +
        $"\t{Og}" +
        $"\t{Ess}" +
        $"\t{MutCount}" +
        $"\t{Hemizygosity}" +
        $"\t{Nullizygosity}" +
        $"\t{Mixture}";

    public static double CalcPloidy(Karyotype kar, GenRef genRef)
        => 2.0 * kar.GenomeLen() / genRef.GetGenomeLen(kar.Sex);
    
    public static SampleStat GetSampleStat(Sample sample, GenRef genRef, FitParams fParams)
    {
        var kar = sample.Karyotype;

        double ploidy = CalcPloidy(kar, genRef);
        double stress = Fitness.StressTerm(genRef.GetGenomeLen(kar.Sex), kar.GenomeLen());
        double tsg = Fitness.TsgOgTerm(genRef, kar.GeneCounts[GeneLT.TSG], kar.Sex, fParams.GeneNormalization);
        double og = Fitness.TsgOgTerm(genRef, kar.GeneCounts[GeneLT.OG], kar.Sex, fParams.GeneNormalization);
        double ess = Fitness.EssTerm(kar.GeneCounts[GeneLT.Ess], fParams.GeneNormalization);
        double fitnessVal = 1 + stress * fParams.Stress + (og - tsg) * fParams.TsgOg + ess * fParams.Essentiality;

        double hemizygosity = Fitness.Zygosity(kar.GeneCounts[GeneLT.Ess], 1);
        double nullizygosity = Fitness.Zygosity(kar.GeneCounts[GeneLT.Ess], 0);

        var res = new SampleStat(sample.SampleId, sample.ParentId, kar.Sex, ploidy, fitnessVal, 
            kar.FitnessVal, stress, tsg, og, ess, 
            sample.Events.Count, hemizygosity, nullizygosity, sample.MixtureString());
        return res;
    }
}