// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public record SimParams(
    int Seed,
    SexEnum Sex,
    int EventCount,
    Distribution Distribution,
    GenomeAssembly Assembly,
    FitnessParams Fitness,
    List<Signature>? Signatures = null,
    MCParams? MCParams = null);