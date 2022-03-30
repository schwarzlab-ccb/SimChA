using CommandLine;

namespace SimChA.IO;

public class CmdOptions
{
    // [Option('p', "popSize", Required = false, Default = 1_000_000, HelpText = "The simulation stops when \"p\" clones are created.")]
    // public int StopCount { get; set; }
    //
    // [Option('c', "cutoff", Required = false, Default = 0, HelpText = "Minimal fraction of a subclone to be considered for an output.")]
    // public float CutOff { get; set; }
    //
    // [Option('r', "reps", Required = false, Default = 1, HelpText = "Number of repetitions")]
    // public int RepsCount { get; set; }
    
    [Option('O', "output", Required = false, Default = "../out", HelpText = "The path to the output files.")]
    public string OutputPath { get; set; }
    
    [Option('C', "config", Required = false, Default = "", HelpText = "A json file with configuration of the experiment.")]
    public string ConfigFile { get; set; }
    
    [Option('N',  Required = false, Default = false, HelpText = "Use newline in logs (useful for batch execution)")]
    public bool Newline { get; set; }
    
    
    // [Option('s', "seed", Required = false, Default = -1, HelpText = "Random seed. If not specified, one is selected at random.")]
    // public int Seed { get; set; }
}