namespace SimChA.DataTypes;

// Parameters used in the Optimization of SimChA Parameters
// Input: 


public record OptimizationParams(
    //int NumSamplesMin, // Burn-in samples
    string Mode, // "Events" or "Fitness"
    int NumSamplesTotal, // Number of true samples
    double StepSize, // Maximum percentage change for a parameter with each iteration
    int MaxTries, // Maximum number of attempts at getting a new event parameter
    bool WriteIntermediate, // Write intermediate better-scoring parameter set to file
    int WriteFrequency, // How often to write intermediate parameter set to file
    double AcceptanceFactor, // Factor for accepting new parameter set
    bool EventWeightsOnly, // Whether to only optimize event weights. If false, optimize event lengths too
    bool ResetSeed, // Whether to reset the random seed for each proposed parameter set
    bool UseABC, // Whether to use Approximate Bayesian Computation
    bool UseSegLength, // Whether to use segment length in optimization distance
    string SegLengthType, // What kind of segment length distribution to use in optimization distance
    // Stratified, All, or Mean
    bool SegmentCountWeighted, // Whether to weight segment count histograms in optimization distance
    bool UsePloidy, // Whether to use ploidy in optimization distance
    bool UseBreakpoints, // Whether to use breakpoints in optimization distance
    bool BreakpointsPerChrom, // If using breakpoints, count per chromosome or per genome
    long BreakpointsBinSize, // Size of bins for breakpoints if BreakpointsPerChrome is false
    int ParamVariationMode, // How many parameters to vary at each step. 
    // 0 = vary all, 1 = vary one, 2 = vary two, etc.
    int GCInterval, // Interval for Garbage Collection (how many simulated runs before performing GC)
    double PloidyCutoff, // Cutoff for ploidy in optimization distance
    int SegLengthCutoff, // Cutoff for segment length in optimization distance
    bool LogTransformSegLength, // Whether to log-transform segment length in optimization distance
    string OptimizationMethod, // Method for optimization - MetropolisHastings, or SimulatedAnnealing, AdaptiveSimulatedAnnealing, or StepSizeDecay
    double StartTemp, // Starting temperature for simulated annealing
    double CoolingRate, // Cooling rate (must be less than 1) used in adapative simulated annealing or alternatively the decay rate for StepSizeDecay
    string StepSizeDecayType, // Type of decay for step size in StepSizeDecay - Linear, Exponential, or Inverse
    double MinStepSize, // Minimum step size for optimization
    bool BoundedLengths, // Whether to bound segment lengths in optimization
    int MaxLength, // Maximum length parameter
    int MinLength, // Minimum length parameter
    bool UseHomozygousDeletion, // Flag to use homozygous deletions in fitness optimization distance
    bool UseCNAlongGenome // Flat to use copy number along the genome in fitness optimization distance
    );