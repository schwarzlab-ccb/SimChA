namespace SimChA.DataTypes;

// Parameters used in the Evolution Modeling
// Input: 

public record EvoParams(
    bool WithFitness, // Evolve with or without fitness
    int KSteps, // Number of steps to evolve between a time period
    Distribution StepDistribution, // Distribution of steps (Exponential or Unit for fixed)
    bool EvolveInTime, // Evolve the clone in time
    double MutationRate, // Mutation rate
    int NumIterations, // Number of true samples
    double ThetaFitness, // exponential multiplier for fitness
    bool SimulatedAnnealing, // Flag to turn on and off simulated annealing of evolution
    double Temperature, // Initial temperature for simulated annealing
    double CoolingRate, // Cooling rate for simulated annealing
    int MaxTries // Maximum number of tries to generate a new set of events in a time step
    );