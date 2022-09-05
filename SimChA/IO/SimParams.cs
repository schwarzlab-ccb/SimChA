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