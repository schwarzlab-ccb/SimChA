// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

[Serializable]
public record CNEventP(CNEventType Type, double Prob, Dictionary<string, double>? Params = null);