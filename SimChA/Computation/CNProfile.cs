// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.Computation;

public abstract class CNProfile
{
    public static double CalcPloidy(Karyotype kar, GenRef genRef)
        => 2.0 * kar.GenomeLen() / genRef.GetGenomeLen(kar.Sex);
    
    public static double CalcCoverage(Karyotype kar, GenRef genRef) 
    =>  (genRef.GetGenomeLen(kar.Sex, false) - kar.MissingLen()) / (double) genRef.GetGenomeLen(kar.Sex,false);
    
    public static CloneStat GetCloneStats(Sample sample, CloneIn clone, GenRef genRef, FitnessParams fParams, Dictionary<string, Karyotype> karMap)
    {
        var kar = karMap[clone.CloneId];

        double ploidy = CalcPloidy(kar, genRef);
        double coverage = CalcCoverage(kar, genRef);
        
        var tsgCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.TumorSuppressor], kar);
        var ogCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.Oncogene], kar);
        var essCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.Essentiality], kar);

        
        double stress = Fitness.StressTerm(genRef.GetGenomeLen(kar.Sex), kar.GenomeLen());
        double tsg = -Fitness.TsgOgTerm(genRef, tsgCNs, kar.Sex, fParams.NormalizeGenes);
        double og = Fitness.TsgOgTerm(genRef, ogCNs, kar.Sex, fParams.NormalizeGenes);
        double ess = Fitness.EssTerm(genRef, essCNs, kar.Sex, fParams.NormalizeGenes, fParams.Haploinsufficiency);
        double fitness = Fitness.CalculateFromComponents(stress, tsg+og, ess, fParams);

        // Get mutation count
        // sample.EventDescs might be empty
        var found_events = sample.EventDescs.TryGetValue(clone.CloneId, out var events);
        int mutCount = found_events ? events.Last().Depth : 0;
        
        double hemizygosity = Fitness.Zygosity(genRef, essCNs, 1);
        double nullizygosity = Fitness.Zygosity(genRef, essCNs, 0);

        return new CloneStat(sample.SampleId, clone.CloneId, ploidy, coverage, fitness, clone.FitnessTarget, stress, tsg, og, ess, mutCount, hemizygosity, nullizygosity);
    }
}
