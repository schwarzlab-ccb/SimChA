using SimChA.IO;

namespace SimChA.Data;

public record SampleStat(
    string SampleId,
    string ParentId,
    SexType Sex,
    double Ploidy,
    double Fitness,
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
        $"\t{Fitness}" +
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

        var tsgCNs = Computation.Fitness.CalcCNs(genRef.GeneLists[GeneListType.TumorSuppressor], kar);
        var ogCNs = Computation.Fitness.CalcCNs(genRef.GeneLists[GeneListType.Oncogene], kar);
        var essCNs = Computation.Fitness.CalcCNs(genRef.GeneLists[GeneListType.Essentiality], kar).ToList();

        double stress = Computation.Fitness.StressTerm(genRef.GetGenomeLen(kar.Sex), kar.GenomeLen());
        double tsg = -Computation.Fitness.TsgOgTerm(genRef, tsgCNs, kar.Sex, fParams.GeneNormalization);
        double og = Computation.Fitness.TsgOgTerm(genRef, ogCNs, kar.Sex, fParams.GeneNormalization);
        double ess = Computation.Fitness.EssTerm(genRef, essCNs, kar.Sex, fParams.GeneNormalization);
        double fitness = Computation.Fitness.CalculateFromComponents(stress, tsg + og, ess, fParams);

        double hemizygosity = Computation.Fitness.Zygosity(genRef, essCNs, 1);
        double nullizygosity = Computation.Fitness.Zygosity(genRef, essCNs, 0);

        var res = new SampleStat(sample.SampleId, sample.ParentId, kar.Sex, ploidy, fitness, 
            kar.FitnessVal, stress, tsg, og, ess, 
            sample.Events.Count, hemizygosity, nullizygosity, sample.MixtureString());
        return res;
    }
}