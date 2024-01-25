namespace SimChA.DataTypes;

// Parameters used in the Optimization of SimChA Parameters
// Input: 


public record OptimizationParams(
    //int NumSamplesMin, // Burn-in samples
    string Mode, // "Events" or "Fitness"
    int NumSamplesTotal, // Number of true samples
    double StepFactor, // Maximum percentage changewith each iteration
    int MaxTries, // Maximum number of attempts at getting a new event parameter
    bool WriteIntermediate, // Write intermediate better-scoring parameter set to file
    int WriteFrequency, // How often to write intermediate parameter set to file
    double AcceptanceFactor, // Factor for accepting new parameter set
    bool BreakpointsPerChrom, // Whether to use breakpoints per chromosome or per genome
    long BreakpointsBinSize, // Size of bins for breakpoints if BreakpointsPerChrome is false
    bool EventWeightsOnly // Whether to only optimize event weights. If false, optimize event lengths too
    );