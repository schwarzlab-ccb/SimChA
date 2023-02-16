// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

namespace SimChA.DataTypes;

public record CNEvent(
    string CloneName, 
    CNEventType AberrationType, 
    int NrOfMutation,
    string Region = "",
    double DeltaFitness = 0, 
    double TotalFitness = 0);