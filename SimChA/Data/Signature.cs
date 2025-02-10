using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;

namespace SimChA.Data;

[Serializable]
public record Signature(double Prob, List<CNEventPars> Events) : IHasProb;

