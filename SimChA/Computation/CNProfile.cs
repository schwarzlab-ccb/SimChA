// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.Computation;

public abstract class CNProfile
{
    public static CloneStat GetCloneStats(Sample sample, CloneIn clone, GenRef genRef, FitnessParams fParams, Dictionary<int, Karyotype> karMap)
    {
        var kar = karMap[clone.CloneId];

        double ploidy = kar.CalcPloidy();
        double coverage = kar.CalcCoverage();
        
        var tsgCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.TumorSuppressor], kar);
        var ogCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.Oncogene], kar);
        var essCNs = Fitness.CalcCNs(genRef.GeneLists[GeneListType.Essentiality], kar);
        
        double fitness = Fitness.Calculate(kar, genRef, fParams);
        double stress = Fitness.StressTerm(kar.GenomeLen(), kar.SexXX);
        double tsg = -Fitness.TsgOgTerm(tsgCNs, kar.SexXX);
        double og = Fitness.TsgOgTerm(ogCNs, kar.SexXX);
        double ess = Fitness.EssTerm(essCNs);

        return new CloneStat(sample.SampleId, clone.CloneId, ploidy, coverage, fitness, stress, tsg, og, ess);
    }
}