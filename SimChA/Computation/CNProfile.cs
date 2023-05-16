// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.Computation;

public abstract class CNProfile
{
    public static CloneStat GetCloneStats(CloneIn clone, Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> geneLists, FitnessParams fParams, Dictionary<int, Karyotype> karMap)
    {
        var kar = karMap[clone.CloneId];
        var tsgCNs = Fitness.CalcCNs(geneLists[GeneListType.TumorSuppressor], kar);
        var ogCNs = Fitness.CalcCNs(geneLists[GeneListType.Oncogene], kar);
        var essCNs = Fitness.CalcCNs(geneLists[GeneListType.Essentiality], kar);
        
        double fitness = Fitness.Calculate(kar, geneLists, fParams);
        double stress = Fitness.StressTerm(kar.GenomeLen(), kar.SexXX);
        double tsg = -Fitness.TsgOgTerm(tsgCNs, kar.SexXX);
        double og = Fitness.TsgOgTerm(ogCNs, kar.SexXX);
        double ess = Fitness.EssTerm(essCNs);

        return new CloneStat(clone.CloneId, fitness, stress, tsg, og, ess);
    }
}