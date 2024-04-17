using CommandLine;
using SimChA.DataTypes;

namespace SimChA.IO;

public class CmdOptions
{
    [Option('O', "output", Required = false, Default = "./out", HelpText = "The path to the output files.")]
    public string OutputPath { get; set; }

    [Option('C', "config", Required = false, Default = "./default_params.json", HelpText = "A json file with configuration of the experiment.")]
    public string ConfigFile { get; set; }

    [Option('T', "tree", Required = false, Default = "", HelpText = "A clone file that describes the tree to be built.")]
    public string CloneTreeFile { get; set; }
    
    [Option('R', "repeats", Required = false, Default = 1, HelpText = "Number of repeats to generate if the distance parameter is used.")]
    public int Repeats { get; set; }
    
    [Option('P', "cnprofiles", Required = false, Default = "", HelpText = "File with CNAs, will cause the program to write a scoring file.")]
    public string CNProfiles { get; set; }

    [Option('D', "data", Required = false, Default = "./data/hg19", HelpText = "The path for the tsv-files with the chromosome list, and the essential-, tsg- and og-genes with scores." +
                                                                               " The files should be named: chromosomes.tsv, essential.tsv, tsg.tsv, og.tsv and contained in a folder with the same name as the assembly used.")]
    public string DataFolder { get; set; }

    [Option('M', "mcmc", Required = false, Default = false, HelpText = "Run the Markov Chain Monte Carlo simulation of mutational events. The argument is a path to a file that lists the fitness of individual clones.")]
    public bool UseMCMC { get; set; }
    
    [Option('s', Required = false, Default = false, HelpText = "Calculate consistent copy numbers segmentation. The output file, consistent_CNs.tsv, will have NA if the original sample did not have data in a given region.")]
    public bool CalcConsistentCNs { get; set; }

    [Option('V', "variants", Required = false, Default = false, HelpText = "Use the included FASTA reference files for SNVs.")]
    public bool UseVariants { get; set; }

    [Option('F', "fasta", Required = false, Default = false, HelpText = "Produce an output FASTA file of the final simulated karyotype, based on the input reference genome.")]
    public bool WriteFasta { get; set; }

    [Option('B', "bootstrap", Required = false, Default = "", HelpText = "Bootstrap sampling of provided fitness from the input file.")]
    public string BootstrapFile { get; set; }

    [Option("optimization", Required = false, Default = false, HelpText = "Run the optimization model.")]
    public bool RunOptimization { get; set; }
    [Option("bin-samples", Required = false, Default = false, HelpText = "Use SimChA to bin copy-number profiles into 1MB-sized bins")]
    public bool BinSamples { get; set; }
    [Option("binned-samples", Required = false, Default = "", HelpText = "Path to the binned data")]
    public string BinnedSamples { get; set; }
    [Option("autosomes-only", Required = false, Default = false, HelpText = "Only consider autosomes for fitness calculations")]
    public bool AutosomesOnly { get; set; }

    [Option("target-params", Required = false, Default = "", HelpText = "A json file with the target set of parameters for parameter inference.")]
    public string TargetParams { get; set; }
    [Option("event-counts", Required = false, Default = "", HelpText = "A tsv file with the event counts for each sample for parameter inference.")]
    public string EventCounts { get; set; }
    public ExecMode ExecMode
    {
        get
        {
            if (BinSamples)
            {
                return ExecMode.BinSamples;
            }
            if (RunOptimization)
            {
                return ExecMode.RunOptimization;
            }
            if (CloneTreeFile != "")
            {
                return ExecMode.Tree;
            }
            if (CNProfiles != "")
            {
                return ExecMode.Profiles;
            }
            if (BootstrapFile != "")
            {
                return ExecMode.Bootstrap;
            }
            if (UseVariants || WriteFasta)
            {
                return ExecMode.ParseGenContents;
            }
            return ExecMode.None;
        }
    }
}