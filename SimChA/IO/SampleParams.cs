using SimChA.Computation;
using SimChA.Data;
using SimChA.DataTypes;

namespace SimChA.IO;

public record SampleParams(
    SexType Sex = SexType.None,
    DistType EventDist = DistType.Uniform,
    double EventMean = 1,
    DistType FitDist = DistType.Uniform,
    double FitMean = 1,
    bool AutosomesOnly = false,
    bool TetraploidStart = false
);