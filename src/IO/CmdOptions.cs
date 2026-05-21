using CommandLine;

namespace SimChA.IO;

public class CmdOptions
{
    [Option('O', "output", Required = false, HelpText = "The path to the output files.")]
    public string OutputPath { get; set; } = "./out";

    [Option('C', "config", Required = false, HelpText = "A json file with configuration of the experiment.")]
    public string ConfigFile { get; set; } = "./configs/main_config.json";

    [Option('T', "tree", Required = false, HelpText = "A clone file that describes the tree to be built.")]
    public string CloneTreeFile { get; set; } = "";

    [Option('R', "repeats", Required = false,
        HelpText = "Number of repeats to generate if the distance parameter is used.")]
    public int Repeats { get; set; } = 1;
    
    [Option('P', "cnprofiles", Required = false, HelpText = "File with CNAs, will cause the program to write a scoring file.")]
    public string CNProfiles { get; set; } = "";
    
    [Option('e', "evolution-mode", Required = false, Default = false, HelpText = "Flag to execute evolution mode. In this mode, events are selected in order to increase fitness.")]
    public bool EvolutionMode { get; set; }

    [Option('m', "match-mode", Required = false, Default = false, HelpText = "Flag to execute fitness matching mode. In this mode, events are selected to minimize the distance to a target fitness.")]
    public bool MatchMode { get; set; }
    
    [Option('s', "segments", Required = false, Default = false, HelpText = "Write out copy numbers segments.")]
    public bool WriteCNs { get; set; }
    
    [Option('S', "consistent-segments", Required = false, Default = false, HelpText = "Write out copy number segments under a minimum consistent segmentation.")]
    public bool WriteConsistentCNs { get; set; }
    

    [Option('v', "variants", Required = false, Default = false, HelpText = "Write out VCF file with the variants of the final simulated karyotype.")]
    public bool WriteVariants { get; set; }

    [Option('f', "fasta", Required = false, Default = false, HelpText = "Write out out a FASTA file for each sample. WARNING! Average file size is 6GB per sample.")]
    public bool WriteFasta { get; set; }

    [Option('z', "zero-index", Required = false, Default = false, HelpText = "Flag for zero-indexed input copy number profiles")]
    public bool ZeroIndexed { get; set; }
    
    [Option("root", Required = false, HelpText = "A path to the folder that will be considered root for relative paths. If not provided, the C# default will be used.")]
    public string RootFolder { get; set; } = ".";

    public ExecMode ExecMode
    {
        get
        {
            if (CloneTreeFile != "" && CNProfiles != "")
            {
                throw new Exception("Cannot run both tree and profiles at the same time.");
            }
            if (CloneTreeFile != "")
            {
                if (Repeats > 1)
                {
                    throw new Exception("Cannot run tree with repeats.");
                }
                return ExecMode.Tree;
            }
            if (CNProfiles != "")
            {
                if (Repeats > 1)
                {
                    throw new Exception("Cannot run profiles with repeats.");
                }
                return ExecMode.Profiles;
            }
            return ExecMode.Repeats;
        }
    }

    public SelectionMode SelectionMode
    {
        get
        {
            if (EvolutionMode && MatchMode)
            {
                throw new Exception("Cannot run both evolution and fitness matching mode at the same time.");
            }
            if (MatchMode)
            {
                return SelectionMode.FitnessMatching;
            }
            if (EvolutionMode)
            {
                return SelectionMode.Evolution;
            }
            return SelectionMode.MonteCarlo;
        }
    }
    
    public bool ShouldParseGenome 
        => WriteVariants || WriteFasta;
    
    public bool Simulate
        => ExecMode != ExecMode.Profiles;
    
    public bool CalcSegments
        => WriteCNs || WriteConsistentCNs;
}
