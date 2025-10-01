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
    double Nullizygosity,
    double Hemizygosity,
    string Mixture,
    int NumRejectedEvents)
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
           "\tnullizygosity" + 
            "\themizygosity" +
           "\tmixture"+
           "\tn_rejected_events";
    
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
        $"\t{Nullizygosity}" +
        $"\t{Hemizygosity}" +
        $"\t{Mixture}"+
        $"\t{NumRejectedEvents}";

    public static double CalcPloidy(Karyotype kar, RefGen refGen)
        => 2.0 * kar.GenomeLen() / refGen.SexGenomeLen[(int) kar.Sex];
    
    public static SampleStat GetSampleStat(Sample sample, RefGen refGen, FitParams fParams)
    {
        var kar = sample.Karyotype;

        double ploidy = CalcPloidy(kar, refGen);
        double stress = Fitness.StressTerm(refGen.SexGenomeLen[(int) kar.Sex], kar.GenomeLen());
        var geneData = refGen.SexGeneLists[(int) kar.Sex];
        double tsg = Fitness.TsgOgTerm(geneData[(int) GeneLT.TSG], kar.GeneCounts[(int) GeneLT.TSG], fParams.GeneNormalization);
        double og = Fitness.TsgOgTerm(geneData[(int) GeneLT.OG], kar.GeneCounts[(int) GeneLT.OG], fParams.GeneNormalization);
        double ess = Fitness.EssTerm(geneData[(int) GeneLT.Ess], kar.GeneCounts[(int) GeneLT.Ess], fParams.GeneNormalization);
        double fitnessVal = 1 + stress * fParams.Stress + (og - tsg) * fParams.TsgOg + ess * fParams.Essentiality;
        
        double nullizygosity = Fitness.Zygosity(geneData[(int) GeneLT.Ess], kar.GeneCounts[(int) GeneLT.Ess], 0);
        double hemizygosity = Fitness.Zygosity(geneData[(int) GeneLT.Ess], kar.GeneCounts[(int) GeneLT.Ess], 1);
        int nRejectedEvents = sample.Events.Sum(e => e.NumRejections);

        var res = new SampleStat(sample.SampleId, sample.ParentId, kar.Sex, ploidy, fitnessVal, 
            kar.FitnessVal, stress, tsg, og, ess, 
            sample.Events.Count, nullizygosity, hemizygosity, sample.MixtureString(), nRejectedEvents);
        return res;
    }
}
