namespace SimChA.DataTypes;

// Parameters used in the Evolution Modeling
// Input: 

public record EvoParams(
    double MutationRate, // Mutation rate
    int NumIterations, // Number of true samples
    double ThetaFitness, // exponential multiplier for fitness
    bool PrintFitnesses// Print the fitnesses of the clone during the MCMC sampling
    );