namespace SimChA.DataTypes;

// Parameters used in the Optimization of SimChA Parameters
// Input: 


public record OptimizationParams(
    //int NumSamplesMin, // Burn-in samples
    string Mode, // "Events" or "Fitness"
    int NumSamplesTotal, // Number of true samples
    double StepFactor, // Maximum percentage changewith each iteration
    int MaxTries, // Maximum number of attempts at getting a new event parameter
    bool WriteIntermediate // Write intermediate better-scoring parameter set to file
    );