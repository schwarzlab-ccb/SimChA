namespace SimChA.IO;

public record SAParams(
    double Acceptance = 0, 
    double EventCost = 0,
    int MaxTries = 1 // Maximum number of tries to generate a new set of events in a time step
);
