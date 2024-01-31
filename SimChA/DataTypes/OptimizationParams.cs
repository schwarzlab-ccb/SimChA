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
    bool EventWeightsOnly, // Whether to only optimize event weights. If false, optimize event lengths too
    bool ResetSeed, // Whether to reset the random seed for each proposed parameter set
    bool UseABC, // Whether to use Approximate Bayesian Computation
    bool UseMeanSeg, // Whether to use mean segment length in optimization distance
    bool UsePloidy, // Whether to use ploidy in optimization distance
    bool UseBreakpoints, // Whether to use breakpoints in optimization distance
    bool BreakpointsPerChrom, // If using breakpoints, count per chromosome or per genome
    long BreakpointsBinSize, // Size of bins for breakpoints if BreakpointsPerChrome is false
    int ParamVariationMode // Mode to vary parameters. 0 = vary all parameters, 
    // 1 = vary a pair of event weights (or a pair of ), 2 = vary only one event at a time
    );