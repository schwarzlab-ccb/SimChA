using SimChA.EventData;

namespace SimChA.DataTypes;

[Serializable]
public record Signature(double Prob, List<CNEventPars> Events) : IHasProb;

