// Created by Felix Schifferdecker, 2022, felix.schifferdecker@gmx.de

using SimChA.DataTypes;

namespace SimChA.IO;

public record Abberation(
    string CloneName, 
    AberrationEnum AberrationType, 
    int NrOfMutation,
    string Region = "",
    double DeltaFitness = 0, 
    double TotalFitness = 0);