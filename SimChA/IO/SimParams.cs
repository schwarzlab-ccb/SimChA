// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.DataTypes;

namespace SimChA.IO;

public record SimParams(
    int Seed,
    bool IsFemale,
    BaseAbbP p_Chromothripsis,
    FractionAbbP p_Inversion,
    FractionAbbP p_Translocation,
    BaseAbbP p_ChromDeletion,
    BaseAbbP p_ChromDuplication,
    FractionAbbP p_InternalDeletion,
    FractionAbbP p_InternalDuplication,
    FractionAbbP p_TailDeletion,
    FractionAbbP p_BreakageFusionBridge,
    BaseAbbP p_WholeGenomeDoubling)
{
    public static SimParams CreateSimParams(int Seed, bool IsFemale, Dictionary<AberrationEnum, BaseAbbP> aberrations)
        => new(
            Seed,
            IsFemale,
            aberrations[AberrationEnum.Chromothripsis],
            (FractionAbbP)aberrations[AberrationEnum.Inversion],
            (FractionAbbP)aberrations[AberrationEnum.Translocation],
            aberrations[AberrationEnum.ChromDeletion],
            aberrations[AberrationEnum.ChromDuplication],
            (FractionAbbP)aberrations[AberrationEnum.InternalDeletion],
            (FractionAbbP)aberrations[AberrationEnum.InternalDuplication],
            (FractionAbbP)aberrations[AberrationEnum.TailDeletion],
            (FractionAbbP)aberrations[AberrationEnum.BreakageFusionBridge],
            aberrations[AberrationEnum.WholeGenomeDoubling]);
}


    // SimParams(int seed)
    // {
    //     Seed = seed;
    //     IsFemale = true;
    //     P_Chromothripsis = new BaseAbbP(1);
    //     P_Inversion = new FractionAbbP( 1, .01);
    //     P_Translocation = new FractionAbbP(1, .1);
    //     P_ChromDeletion = new BaseAbbP( 1);
    //     P_ChromDuplication = new BaseAbbP(1);
    //     P_InternalDeletion = new FractionAbbP(1, .01);
    //     P_InternalDuplication = new FractionAbbP( 1, .01);
    //     P_TailDeletion = new FractionAbbP(1, .1);
    //     P_BreakageFusionBridge = new FractionAbbP(1, .1);
    //     P_WholeGenomeDoubling = new BaseAbbP( 1);
    // }
    //
    // public List<BaseAbbP> AberrationsList() => new List<BaseAbbP>
    // {
    //     P_Chromothripsis, P_Inversion, P_Translocation, P_ChromDeletion, P_ChromDuplication, P_InternalDeletion, P_InternalDuplication, P
    // }
    //
    // public double SumRates() => AberrationsList().Sum(ar => ((BaseAbbP) ar.Value).Likelihood);
    //
    // public static Dictionary<AberrationEnum, object> DefaultAberrations()
    //     => new()
    //     {
    //         [AberrationEnum.Chromothripsis] = new BaseAbbP(1),
    //         [AberrationEnum.Inversion] = new FractionAbbP( 1, .01),
    //         [AberrationEnum.Translocation] = new FractionAbbP(1, .1),
    //         [AberrationEnum.ChromDeletion] = new BaseAbbP( 1),
    //         [AberrationEnum.ChromDuplication] = new BaseAbbP(1),
    //         [AberrationEnum.InternalDeletion] = new FractionAbbP(1, .01),
    //         [AberrationEnum.InternalDuplication] = new FractionAbbP( 1, .01),
    //         [AberrationEnum.TailDeletion] = new FractionAbbP(1, .1),
    //         [AberrationEnum.BreakageFusionBridge] = new FractionAbbP(1, .1),
    //         [AberrationEnum.WholeGenomeDoubling] = new BaseAbbP( 1),
    //     };

