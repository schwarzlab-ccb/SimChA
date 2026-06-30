// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using System.Reflection;
using SimChA.Data;

namespace SimChA.IO;

public record SimChAConfig(
    SimParams SimParams,
    FitParams FitParams,
    List<Signature>? Signatures = null,
    EvoParams? EvoParams = null,
    // Top-level root folder used to resolve relative paths. Read from the config if present and
    // overridable on the command line; stamped with the effective value in the output config.
    string? Root = null,
    // The SimChA version that produced this config. Set on output; ignored on input.
    string? Version = null)
{
    // Reads the config file referenced by the command-line options and resolves the effective root,
    // assembly and gene set with precedence command line > config > default. The root is applied as
    // the working directory so relative paths resolve against it, and the effective values plus the
    // running version are stamped into the returned config so the output records exactly what was used.
    public static SimChAConfig Load(CmdOptions options)
    {
        // The command-line root is applied first so the config path resolves against it.
        string? cmdRoot = string.IsNullOrEmpty(options.RootFolder) ? null : Path.GetFullPath(options.RootFolder);
        if (cmdRoot != null)
        {
            Environment.CurrentDirectory = cmdRoot;
        }

        var config = FileIO.ReadSimChAConfig(options.ConfigFile);

        string effectiveRoot = cmdRoot
            ?? (string.IsNullOrEmpty(config.Root) ? Path.GetFullPath(".") : Path.GetFullPath(config.Root));
        Environment.CurrentDirectory = effectiveRoot;

        // Effective assembly / gene set with precedence: command line > config > config default.
        string effectiveAssembly = string.IsNullOrEmpty(options.AssemblyFolder)
            ? config.SimParams.Assembly
            : options.AssemblyFolder;
        string effectiveGeneSet = string.IsNullOrEmpty(options.GeneSetFolder)
            ? config.FitParams.GeneSet
            : options.GeneSetFolder;
        string version = Assembly.GetExecutingAssembly().GetName().Version?.ToString() ?? "unknown";

        return config with
        {
            SimParams = config.SimParams with { Assembly = effectiveAssembly },
            FitParams = config.FitParams with { GeneSet = effectiveGeneSet },
            Root = effectiveRoot,
            Version = version
        };
    }
}
