// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.Computation;

public abstract class CNProfile
{
    public static double CalcPloidy(Karyotype kar, GenRef genRef)
        => 2.0 * kar.GenomeLen() / genRef.GetGenomeLen(kar.SexXX);
    
    public static double CalcCoverage(Karyotype kar, GenRef genRef) 
    =>  (genRef.GetGenomeLen(kar.SexXX) - kar.MissingLen()) / (double) genRef.GetGenomeLen(kar.SexXX);
    
    public static CloneStat GetCloneStats(Sample sample, CloneIn clone, GenRef genRef, FitnessParams fParams, Dictionary<int, Karyotype> karMap)
    {
        var kar = karMap[clone.CloneId];

        double ploidy = CalcPloidy(kar, genRef);
        double coverage = CalcCoverage(kar, genRef);
        
        var tsgCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.TumorSuppressor], kar);
        var ogCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.Oncogene], kar);
        var essCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.Essentiality], kar);
        
        double fitness = Fitness.Calculate(kar, genRef, fParams);
        double stress = Fitness.StressTerm(genRef.GetGenomeLen(kar.SexXX), kar.GenomeLen());
        double tsg = -Fitness.TsgOgTerm(genRef, tsgCNs, kar.SexXX);
        double og = Fitness.TsgOgTerm(genRef, ogCNs, kar.SexXX);
        double ess = Fitness.EssTerm(genRef, essCNs, kar.SexXX);

        return new CloneStat(sample.SampleId, clone.CloneId, ploidy, coverage, fitness, stress, tsg, og, ess);
    }
}
