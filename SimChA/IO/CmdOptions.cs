using CommandLine;

namespace SimChA.IO;

public class CmdOptions
{
    [Option('s', "stopcount", Required = false, Default = 1000, HelpText = "The simulation stops when \"s\" clones are created.")]
    public int StopCount { get; set; }
    
    [Option('c', "cutoff", Required = false, Default = 1, HelpText = "The minimum number of cells to be written to output.")]
    public int CutOff { get; set; }
    
    [Option('O', "output", Required = false, Default = "out", HelpText = "The path to the output files.")]
    public string OutputPath { get; set; }
    
    [Option('D', "dotGraph", Required = false, Default = "SimChA.dot", HelpText = "The path to the DOT graph output file.")]
    public string DotPath { get; set; }
}