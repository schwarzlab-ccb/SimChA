using SimChA.Computation;
using SimChA.Data;

namespace SimChA.IO;

public record SimParams(
    int Seed = 0,
    SexType Sex = SexType.Any,
    DistType RateDist = DistType.Uniform,
    double RateMean = 1,
    DistType FitDist = DistType.Uniform,
    double FitMean = 1,
    bool TetraploidStart = false,
    MixtureType Mixture = MixtureType.Constant
);