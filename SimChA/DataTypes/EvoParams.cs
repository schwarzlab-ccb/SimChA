namespace SimChA.DataTypes;

// Parameters used in the Evolution Modeling

public record EvoParams(
    bool WithFitness, // Evolve with or without fitness
    int MaxTries, // Maximum number of tries to generate a new set of events in a time step
    double EventCost, // Flag to turn on and off event parameters being weighted by ploidy
    bool TetraploidStart, // Flag to start the simulation from a tetraploid state
    bool EvolveInTime, // Flag to evolve in time
    float MaxTime, // Maximum time to evolve
    float MutationRate // Average number of mutations per unit time
);
