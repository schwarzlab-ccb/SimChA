// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.IO;

namespace SimChA.Simulation;

public class AberrationsInfo
{
    public Dictionary<AberrationEnum, BaseAbbP> Map { get; }
    public double RatesSum { get; }

    public AberrationsInfo(SimParams simParams)
    {
        Map = new Dictionary<AberrationEnum, BaseAbbP>
        {
            [AberrationEnum.Chromothripsis] = simParams.p_Chromothripsis,
            [AberrationEnum.Inversion] = simParams.p_Inversion,
            [AberrationEnum.Translocation] = simParams.p_Translocation,
            [AberrationEnum.ChromDeletion] = simParams.p_ChromDeletion,
            [AberrationEnum.ChromDuplication] = simParams.p_ChromDuplication,
            [AberrationEnum.InternalDeletion] = simParams.p_InternalDuplication,
            [AberrationEnum.InternalDuplication] = simParams.p_InternalDuplication,
            [AberrationEnum.TailDeletion] = simParams.p_TailDeletion,
            [AberrationEnum.BreakageFusionBridge] = simParams.p_BreakageFusionBridge,
            [AberrationEnum.WholeGenomeDoubling] = simParams.p_WholeGenomeDoubling
        };
        RatesSum = Map.Sum(ar => ar.Value.Likelihood);
    }

    public static Dictionary<AberrationEnum, BaseAbbP> DefaultAberrations()
        => new()
        {
            [AberrationEnum.Chromothripsis] = new BaseAbbP(1),
            [AberrationEnum.Inversion] = new FractionAbbP( 1, .01),
            [AberrationEnum.Translocation] = new FractionAbbP(1, .1),
            [AberrationEnum.ChromDeletion] = new BaseAbbP( 1),
            [AberrationEnum.ChromDuplication] = new BaseAbbP(1),
            [AberrationEnum.InternalDeletion] = new FractionAbbP(1, .01),
            [AberrationEnum.InternalDuplication] = new FractionAbbP( 1, .01),
            [AberrationEnum.TailDeletion] = new FractionAbbP(1, .1),
            [AberrationEnum.BreakageFusionBridge] = new FractionAbbP(1, .1),
            [AberrationEnum.WholeGenomeDoubling] = new BaseAbbP( 1),
        };
}