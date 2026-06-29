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

    [Option("seed-karyotypes", Required = false,
        HelpText = "File with sample_id and karyotype columns used to seed the simulation from exact karyotype structure.")]
    public string SeedKaryotypes { get; set; } = "";
    
    [Option('m', "mode", Required = false, Default = "evolution",
        HelpText = "The event selection mode: 'basic' (events are selected at random), " +
                   "'evolution' (events are selected to increase fitness), or " +
                   "'matching' (events are selected to minimize the distance to a target fitness).")]
    public string Mode { get; set; } = "evolution";
    
    [Option('s', "segments", Required = false, Default = false, HelpText = "Write out copy numbers segments.")]
    public bool WriteCNs { get; set; }

    [Option('p', "profile", Required = false, Default = false, HelpText = "Start with given copy number profile.")]
    public bool StartWithProfile { get; set; }
    
    [Option('S', "consistent-segments", Required = false, Default = false, HelpText = "Write out copy number segments under a minimum consistent segmentation.")]
    public bool WriteConsistentCNs { get; set; }

    [Option('K', "karyotypes", Required = false, Default = false, HelpText = "Write out karyotypes of each sample.")]
    public bool WriteKaryotypes { get; set; }

    [Option('v', "variants", Required = false, Default = false, HelpText = "Write out VCF file with the variants of the final simulated karyotype. Requires data/<assembly>/genome.fa (see scripts/DownloadRefData.sh).")]
    public bool WriteVariants { get; set; }

    [Option('f', "fasta", Required = false, Default = false, HelpText = "Write out a FASTA file for each sample. Requires data/<assembly>/genome.fa (see scripts/DownloadRefData.sh). WARNING! Average file size is 6GB per sample.")]
    public bool WriteFasta { get; set; }

    [Option('z', "zero-index", Required = false, Default = false, HelpText = "Flag for zero-indexed input copy number profiles")]
    public bool ZeroIndexed { get; set; }
    
    [Option('r', "root", Required = false, HelpText = "A path to the folder that will be considered root for relative paths. If not provided, the C# default will be used.")]
    public string RootFolder { get; set; } = ".";
    
    [Option('d', "delta", Required = false, Default = false, HelpText = "Will also print the changes caused by events.")]
    public bool Debug { get; set; }
    
    [Option('k', "karyotype",  Required = false, Default = false, HelpText = "Will also print the karyotype after each event.")]
    public bool Karyotype { get; set; }

    // True when a CN profile should seed the simulation (-P together with -p)
    // rather than being scored only.
    public bool Seeded => StartWithProfile && CNProfiles != "";

    public bool SeededFromKaryotype => SeedKaryotypes != "";

    public ExecMode ExecMode
    {
        get
        {
            if (Seeded && SeededFromKaryotype)
            {
                throw new Exception("Cannot seed from both a CN profile and a karyotype file at the same time.");
            }
            if (CloneTreeFile != "" && CNProfiles != "" && !Seeded)
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
            if (CNProfiles != "" && !Seeded)
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
        => Mode.ToLowerInvariant() switch
        {
            "basic" => SelectionMode.MonteCarlo,
            "evolution" => SelectionMode.Evolution,
            "matching" => SelectionMode.FitnessMatching,
            _ => throw new Exception(
                $"Unknown mode '{Mode}'. Valid values are: basic, evolution, matching.")
        };
    
    public bool ShouldParseGenome 
        => WriteVariants || WriteFasta;
    
    public bool Simulate
        => ExecMode != ExecMode.Profiles;
    
    public bool CalcSegments
        => WriteCNs || WriteConsistentCNs;
}
