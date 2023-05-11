// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

namespace SimChA.EventData;

public record CNEventDesc(
    int CloneId,  
    CNEventType EventType,
    string Description = "",
    double DeltaFitness = 0, 
    double TotalFitness = 0);