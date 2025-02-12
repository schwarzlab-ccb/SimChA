using SimChA.Computation;
using SimChA.DataTypes;

namespace SimChA.EventData;

[Serializable]
public record CNEventPars(CNEventType Type, double Prob, long Size = 0, double Frag = 0) : IHasProb;
