using SimChA.Data;
using SimChA.IO;

namespace SimChA.Computation;

public abstract class CNProfile
{
    public static double CalcPloidy(Karyotype kar, GenRef genRef)
        => 2.0 * kar.GenomeLen() / genRef.GetGenomeLen(kar.Sex);
    
    public static double CalcCoverage(Karyotype kar, GenRef genRef) 
    =>  (genRef.GetGenomeLen(kar.Sex, false) - kar.MissingLen()) / (double) genRef.GetGenomeLen(kar.Sex,false);
    
    public static CloneStat GetCloneStats(string sampleId, Sample sample, GenRef genRef, FitParams fParams)
    {
        var kar= sample.Karyotype;

        double ploidy = CalcPloidy(kar, genRef);
        double coverage = CalcCoverage(kar, genRef);
        
        var tsgCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.TumorSuppressor], kar);
        var ogCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.Oncogene], kar);
        var essCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.Essentiality], kar);
        
        double stress = Fitness.StressTerm(genRef.GetGenomeLen(kar.Sex), kar.GenomeLen());
        double tsg = -Fitness.TsgOgTerm(genRef, tsgCNs, kar.Sex, fParams.GeneNormalization);
        double og = Fitness.TsgOgTerm(genRef, ogCNs, kar.Sex, fParams.GeneNormalization);
        double ess = Fitness.EssTerm(genRef, essCNs, kar.Sex, fParams.GeneNormalization);
        double fitness = Fitness.CalculateFromComponents(stress, tsg+og, ess, fParams);
        
        double hemizygosity = Fitness.Zygosity(genRef, essCNs, 1);
        double nullizygosity = Fitness.Zygosity(genRef, essCNs, 0);

        var res = new CloneStat(sampleId, sample.SampleId, ploidy, coverage, fitness, kar.FitnessVal, stress, tsg, og, ess, sample.Events.Count, hemizygosity, nullizygosity);
        return res;
    }
}
