namespace SimChA.DataTypes;

// Parameters used in the Evolution Modeling
// Input: 

public record EvoParams(
    bool WithFitness, // Evolve with or without fitness
    bool EventBlock, // Flag on whether to sample events in blocks (Poisson/geometric distributed)
    Distribution StepDistribution, // Distribution of steps (Exponential or Unit for fixed)
    bool EvolveInTime, // Evolve the clone in time
    double MutationRate, // Mutation rate
    double ThetaFitness, // exponential multiplier for fitness
    int MaxTries, // Maximum number of tries to generate a new set of events in a time step
    bool EventCost, // Flag to turn on and off event parameters being weighted by ploidy
    bool DynamicMutRate, // Flag to turn on and off dynamic mutation rate weighted by ploidy
    double MaxTime, // Maximum time to evolve each clone
    bool ContinuousTime, // Flag to turn on and off continuous time evolution
    bool TetraploidStart, // Flag to start the simulation from a tetraploid state
    double WGDAccelerationFactor, // Factor to scale the mutation rates of WGD+ samples
    double ChromLossEnhancementFactor // Factor to scale the mutation rates of chromosomal loss events, post-WGD
    );
