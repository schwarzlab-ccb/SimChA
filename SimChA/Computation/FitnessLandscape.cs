using SimChA.DataTypes;
using SimChA.IO;
namespace SimChA.Computation;

public static class FitnessLandscape
{
    public static void GenerateFitnessLandscape(GenRef genRef, SimParams simParams, List<Sample> samples, FileIO files)
    {
        GetSampleStats(genRef, simParams, samples);
        int nSteps = 1000;
        var range = Enumerable.Range(0, nSteps).Select(i => i * 1.0 / (nSteps - 1) ).ToList();
        foreach (var sample in samples)
        {
            foreach (var clone in sample.Clones)
            {
                var stats = sample.Stats[clone.CloneId];
                var stress = stats.Stress;
                var tsgog = stats.Tsg + stats.Og;
                var ess = stats.Ess;

                var output = new List<List<double>>();
                foreach (var alpha in range)
                {
                    foreach (var beta in range)
                    {
                        if (alpha + beta > 1)
                        {
                            continue;
                        }
                        var fitness = 1.0 + alpha * stress + beta * tsgog + (1.0 - alpha - beta) * ess;
                        output.Add(new List<double> {alpha, beta, fitness});
                    }
                }
                var filename = $"FitnessLandscape_{sample.SampleId}_{clone.CloneId}.tsv";
                files.WriteFitnessLandscape(filename, output);
            }
        }
    }

    private static void GetSampleStats(GenRef genRef, SimParams simParams, List<Sample> samples)
    {
        foreach (var sample in samples)
        {
            int counter = 1;
            int total = sample.Clones.Count;
            foreach (var clone in sample.Clones)
            {
                Console.Write($"\rSample {sample.SampleId}. Clone {counter++}/{total}.".PadRight(80));
                sample.Stats[clone.CloneId] = CNProfile.GetCloneStats(sample, clone, genRef, simParams.Fitness, sample.Kars);
            }
        }
    }
}