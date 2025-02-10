// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.Data;

namespace SimChA.IO;

public record SimParams(
    int Seed,
    SampleParams SampleParams,
    FitParams FitParams,
    Dictionary<string, Signature>? Signatures = null,
    MHParams? MHParams = null,
    EvoParams? EvoParams = null);