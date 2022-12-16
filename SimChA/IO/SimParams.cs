// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.DataTypes;

namespace SimChA.IO;

public record SimParams(
    int Seed,
    bool IsFemale,
    float StressFraction,
    float TsgOgFraction,
    float EssentialFraction,
    BaseAbbP p_ChromDeletion,
    BaseAbbP p_ChromDuplication,
    FractionAbbP p_TailDeletion,
    FractionAbbP p_InternalDeletion,
    FractionAbbP p_InternalDuplication,
    FractionAbbP p_InternalInversion,
    BaseAbbP p_Translocation,
    FractionAbbP p_BreakageFusionBridge,
    BaseAbbP p_WholeGenomeDoubling,
    BaseAbbP p_Chromothripsis)
{
    public static SimParams CreateSimParams(int seed, bool isFemale, float stressFraction, float tsgOgFraction, 
        float essentialFraction, Dictionary<AberrationEnum, BaseAbbP> aberrations)
        => new(
            seed,
            isFemale,
            stressFraction,
            tsgOgFraction,
            essentialFraction,
            (BaseAbbP) aberrations[AberrationEnum.ChromDeletion],
            (BaseAbbP) aberrations[AberrationEnum.ChromDuplication],
            (FractionAbbP) aberrations[AberrationEnum.TailDeletion],
            (FractionAbbP) aberrations[AberrationEnum.InternalDeletion],
            (FractionAbbP) aberrations[AberrationEnum.InternalDuplication],
            (FractionAbbP) aberrations[AberrationEnum.InternalInversion],
            (BaseAbbP) aberrations[AberrationEnum.Translocation],
            (FractionAbbP) aberrations[AberrationEnum.BreakageFusionBridge],
            (BaseAbbP) aberrations[AberrationEnum.WholeGenomeDoubling],
            (BaseAbbP) aberrations[AberrationEnum.Chromothripsis]);
}