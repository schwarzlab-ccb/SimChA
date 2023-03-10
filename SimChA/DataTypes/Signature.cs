namespace SimChA.DataTypes;

[Serializable]
public record Signature(string Id, double Prob, double Circularization, List<CNEventP>? Events);