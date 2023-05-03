// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.Computation;

public class CNProfile
{
    public static ProfileStats GetProfileStats(Clone clone, Dictionary<GeneListType, Dictionary<ChrNo, List<Gene>>> geneLists, FitnessParams fParams)
    {
        var kar = clone.Karyotype;
        var tsgCNs = Fitness.CalcCNs(geneLists[GeneListType.TumorSuppressor], kar);
        var ogCNs = Fitness.CalcCNs(geneLists[GeneListType.Oncogene], kar);
        var essCNs = Fitness.CalcCNs(geneLists[GeneListType.Essentiality], kar);
        
        double fitness = Fitness.Calculate(kar, geneLists, fParams);
        double stress = Fitness.StressTerm(kar.GenomeLen(), kar.SexXX);
        double tsg = -Fitness.TsgOgTerm(tsgCNs);
        double og = Fitness.TsgOgTerm(ogCNs);
        double ess = Fitness.EssTerm(essCNs);

        return new ProfileStats(clone.Name, fitness, stress, tsg, og, ess);
    }
}