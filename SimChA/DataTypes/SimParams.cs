// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public record SimParams(
    int Seed,
    bool SexXX,
    GenomeAssembly Assembly,
    FitnessParams Fitness,
    List<Signature>? Signatures);