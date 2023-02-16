// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public record SimParams(
    int Seed,
    bool IsFemale,
    // Multiplicative factors for calculation of fitness 
    float StressFraction,
    float TsgOgFraction,
    float EssentialFraction,
    List<Signature>? Signatures = null);