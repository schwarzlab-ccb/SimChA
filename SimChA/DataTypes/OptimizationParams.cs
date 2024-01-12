namespace SimChA.DataTypes;

// Parameters used in the Optimization of SimChA Parameters
// Input: 


public record OptimizationParams(
    //int NumSamplesMin, // Burn-in samples
    string Mode, // "Events" or "Fitness"
    int NumSamplesTotal, // Number of true samples
    double StepFactor // Maximum percentage changewith each iteration
    );