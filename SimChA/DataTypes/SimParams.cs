// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public record SimParams(
    int Seed,
    SexEnum Sex,
    bool AutosomesOnly,
    double EventCountMean,
    Distribution EventDist,
    FitnessParams Fitness,
    Dictionary<string, Signature>? Signatures = null,
    MCParams? MCParams = null,
    MCTarget? MCTarget = null,
    EvoParams? EvoParams = null);