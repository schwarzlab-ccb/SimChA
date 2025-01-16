// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

namespace SimChA.EventData;

public record CNEventDesc(
    CNEventType EventType,
    int Depth,
    string Description = "",
    double DeltaFitness = 0, 
    double TotalFitness = 0,
    double Time = 0);