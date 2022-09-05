using CommandLine;

namespace SimChA.IO;

public class CmdOptions
{
    [Option('O', "output", Required = false, Default = "../out", HelpText = "The path to the output files.")]
    public string OutputPath { get; set; }

    [Option('C', "config", Required = false, Default = "", HelpText = "A json file with configuration of the experiment.")]
    public string ConfigFile { get; set; }

    [Option('I', Required = true, Default = "", HelpText = "Newick file that describes the tree to be built.")]
    public string NewickFile { get; set; }
}