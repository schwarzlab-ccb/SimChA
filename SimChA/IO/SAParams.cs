namespace SimChA.IO;

public record SAParams(
    double EventCost, // Flag to turn on and off event parameters being weighted by ploidy
    double Acceptance, 
    double EventRate,
    int MaxTries, // Maximum number of tries to generate a new set of events in a time step
    bool EvolveInTime
);
