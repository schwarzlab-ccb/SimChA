// Created by Felix Schifferdecker, 2022, felix.schifferdecker@gmx.de

namespace SimChA.IO;

public record Abberation(
    string CloneName, 
    string AbberationEnum, 
    int NrOfMutation,
    string Region = "",
    float DeltaFitness = 0, 
    float TotalFitness = 0);