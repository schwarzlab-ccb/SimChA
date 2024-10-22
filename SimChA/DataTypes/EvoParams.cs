namespace SimChA.DataTypes;

// Parameters used in the Evolution Modeling
// Input: 

public record EvoParams(
    bool EvolveInTime, // Evolve the clone in time
    double MutationRate, // Mutation rate
    int NumIterations, // Number of true samples
    double ThetaFitness, // exponential multiplier for fitness
    bool SimulatedAnnealing, // Flag to turn on and off simulated annealing of evolution
    double Temperature, // Initial temperature for simulated annealing
    double CoolingRate // Cooling rate for simulated annealing
    );