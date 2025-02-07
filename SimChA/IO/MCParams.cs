namespace SimChA.DataTypes;

// Parameters used in the MCMC sampling
// Input: 

public record MCParams(
    int NumSamplesMin, // Burn-in samples
    int NumSamplesTotal, // Number of true samples
    double ThetaFitness, // exponential multiplier for fitness
    bool MatchFitness, // Match the fitness of the clone to a target fitness
    double SwapEventP, // Probability of completely swapping an event
    double ThresholdFit // Percentage difference allowed between accepted and target fitness
);