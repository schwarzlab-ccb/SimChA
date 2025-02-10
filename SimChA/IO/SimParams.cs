using SimChA.Computation;
using SimChA.Data;

namespace SimChA.IO;

public record SimParams(
    SexType Sex = SexType.Any,
    DistType RateDist = DistType.Uniform,
    double RateMean = 1,
    DistType FitDist = DistType.Uniform,
    double FitMean = 1,
    bool AutosomesOnly = false,
    bool TetraploidStart = false
);