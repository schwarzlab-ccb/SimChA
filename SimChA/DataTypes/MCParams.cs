namespace SimChA.DataTypes;

// Parameters used in the MCMC sampling
// Input: 


public record MCParams(
    int NumBurnIn, int NumSamples, // number of burn-in events, and number of true samples
    double ThetaFitness, // exponential multiplier for fitness
    double SwapEventP // Probability of completely swapping an event
    );