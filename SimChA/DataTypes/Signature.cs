namespace SimChA.DataTypes;

[Serializable]
public record Signature(string Id, double Prob, List<CNEventP> Events);