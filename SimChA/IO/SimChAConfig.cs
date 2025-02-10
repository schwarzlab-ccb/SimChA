// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.Data;

namespace SimChA.IO;

public record SimChAConfig(
    int Seed,
    SimParams SimParams,
    FitParams FitParams,
    bool AutosomesOnly = true,
    List<Signature>? Signatures = null,
    MHParams? MHParams = null,
    EvoParams? EvoParams = null);