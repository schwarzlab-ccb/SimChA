using SimChA.EventData;

namespace SimChA.DataTypes;

[Serializable]
public record Signature(string Id, double Prob, List<CNEventPars> Events);

