namespace SimChA.IO;

public record SAParams(
    double EventCost = 0, // Flag to turn on and off event parameters being weighted by ploidy
    double Acceptance = 0, 
    double EventRate = 1,
    int MaxTries = 1, // Maximum number of tries to generate a new set of events in a time step
    bool EvolveInTime = true
);
