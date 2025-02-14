namespace SimChA.IO;

// Parameters used in the MCMC sampling
// Input: 

public record MHParams(
    int NumSamplesMin = 0, // Burn-in samples
    int NumSamplesTotal = 1, // Number of true samples
    double ThetaFitness = 1, // exponential multiplier for fitness
    bool MatchFitness = true, // Match the fitness of the clone to a target fitness
    double SwapEventP = 1, // Probability of completely swapping an event
    double ThresholdFit = 0 // Percentage difference allowed between accepted and target fitness
);