// Created by Felix Schifferdecker, 2022, felix.schifferdecker@gmx.de

namespace SimChA.DataTypes;

public record Gene(string Name, GenRange Region, double DeltaFitness);