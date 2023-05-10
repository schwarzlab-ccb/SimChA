using CommandLine;

namespace SimChA.IO;

public class CmdOptions
{
    [Option('O', "output", Required = false, Default = "./out", HelpText = "The path to the output files.")]
    public string OutputPath { get; set; }

    [Option('C', "config", Required = false, Default = "./default_params.json", HelpText = "A json file with configuration of the experiment.")]
    public string ConfigFile { get; set; }

    [Option('T', "tree", Required = false, Default = "", HelpText = "Newick file that describes the tree to be built.")]
    public string NewickFile { get; set; }
    
    [Option('R', "repeats", Required = false, Default = 1, HelpText = "Number of repeats to generate if the distance parameter is used.")]
    public int Repeats { get; set; }
    
    [Option('P', "cnprofiles", Required = false, Default = "", HelpText = "File with CNAs, will cause the program to write a scoring file.")]
    public string CNProfiles { get; set; }

    [Option("data", Required = false, Default = "./data", HelpText = "Folder with three files for OGs, TSGs and essential genes in the format of essentials.tsv, ogs.tsv and tsgs.tsv")]
    public string GenesFolder { get; set;}

    [Option('M', "mcmc", Required = false, Default = false, HelpText = "Run the Markov Chain Monte Carlo simulation of mutational events. The argument is a path to a file that lists the fitness of individual clones.")]
    public bool UseMCMC {get; set;}

}