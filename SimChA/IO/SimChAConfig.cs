// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.Data;

namespace SimChA.IO;

public record SimChAConfig(
    SimParams SimParams,
    ChAParams ChAParams,
    FitParams FitParams,
    List<Signature>? Signatures = null,
    MHParams? MHParams = null,
    SAParams? EvoParams = null);