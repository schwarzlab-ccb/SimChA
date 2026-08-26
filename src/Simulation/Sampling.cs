using SimChA.Computation;
using SimChA.EventData;
using MathNet.Numerics.Distributions;
using SimChA.Data;

namespace SimChA.Simulation;

public static class Sampling
{
    public const double FixedBetaAlpha = 0.6;

    public static long GetNormSeg(Random rnd, long contigLen, double meanFrac) 
        => (long) Math.Clamp(Math.Round(contigLen * Normal.Sample(rnd, meanFrac, meanFrac / 3)), 1, contigLen);
    
    public static long GetExpSeg(Random rnd, long contigLen, long meanLen) 
        => (long) Math.Clamp(Math.Round(meanLen * Exponential.Sample(rnd, 1)), 1, contigLen);
    
    public static long GetExpSeg(Random rnd, long contigLen, double meanFrac)
        => (long) Math.Clamp(Math.Round(contigLen * meanFrac * Exponential.Sample(rnd, 1)), 1, contigLen);

    public static long GetParetoSeg(Random rnd, long contigLen, double meanFrac)
        => (long) Math.Clamp(Math.Round(contigLen * SampleParetoLim(rnd, meanFrac)), 1, contigLen);

    public static long GetBetaSeg(Random rnd, long armLength, double meanFrac)
    {
        if (meanFrac >= 1)
        {
            return armLength;
        }
        if (meanFrac <= 0)
        {
            return 1;
        }
        return (long) Math.Clamp(
            Math.Round(armLength * SampleBeta(rnd, meanFrac)),
            1,
            armLength);
    }

    public static long GetPos(Random rnd, long contigLen)
        => rnd.NextInt64(0, contigLen);
    
    public static int GetFragCount(Random rnd, double mean)
        => (int) Math.Max(1, Math.Round(Normal.Sample(rnd, mean, 1)));

    // https://ashpublications.org/blood/article/134/Supplement_1/3767/424006/Chromoplexy-and-Chromothripsis-Are-Important
    public static int GetChromoplexySiteCount(Random rnd)
        => rnd.NextSingle() switch
        {
            var n when n < .46 => 3,
            var n when n < .64 => 4,
            var n when n < .74 => 5,
            var n when n < .79 => 6,
            _ => 2
        };

    public static List<long> GetStopsForShards(Random rnd, long contigLen, int shardCount)
    {
        var stops = new List<long>();
        if (shardCount <= 1)
        {
            return stops;
        }
        for (int i = 1; i < shardCount; i++)
        {
            long newStop = GetPos(rnd, contigLen);
            stops.Add(newStop);
            stops.Sort();
        }
        return stops;
    }
    
    public static IList<double> CreateRandomMixture(Random rnd,  IList<double> concentrations)
        => concentrations.Count != 0 ? new Dirichlet(concentrations.ToArray(), rnd).Sample() : [];

    
    public static IList<double> ConcentrationsToProbabilities(Random rnd, List<double> concentrations, MixtureType mixType)
    {
        return mixType switch
        {
            MixtureType.Single => SelectFromMixture(rnd, concentrations),
            MixtureType.Constant => concentrations,
            MixtureType.Dirichlet => CreateRandomMixture(rnd, concentrations),
            _ => throw new ArgumentOutOfRangeException(nameof(mixType), mixType, null)
        };
    }
    
    public static List<double> SelectFromMixture(Random rnd, IList<double> concentrations)
    {
        int selected = rnd.PickRndIndex(concentrations);
        double[] mixture = new double[concentrations.Count];
        mixture[selected] = 1;
        return mixture.ToList();
    }
    
    public static double SampleContDist(Random rnd, DistType dist, double mean = 1)
    {
        return dist switch
        {
            DistType.Exponential => Exponential.Sample(rnd, 1 / mean),
            DistType.Normal => Normal.Sample(rnd, mean, 0.5),
            _ => mean
        };
    }

    // Gamma is parameterised by (shape, mean) rather than (shape, scale): the shape is a pure
    // dispersion term (CV = 1/sqrt(shape)) so it can be held constant while the mean varies per
    // cancer type, and the scale follows as mean/shape. MathNet takes a rate, i.e. shape/mean.
    // The continuous draw is rounded and clamped to >= 1 so the result is a usable event count,
    // matching Geometric's {1,2,...} support; at the means used here (tens of events) the
    // discretisation shifts the mean by well under a percent.
    public static int SampleDiscDist(Random rnd, DistType dist, double mean, double shape = 1)
    {
        return dist switch
        {
            DistType.Geometric => Geometric.Sample(rnd, 1 / mean),
            DistType.Poisson => Poisson.Sample(rnd, mean),
            DistType.Gamma => Gamma.Sample(rnd, shape, shape / mean).ToCount(),
              _ => (int) mean
        };
    }

    private static int ToCount(this double sample)
        => (int) Math.Max(1, Math.Round(sample));

    private static double GetParetoScale(double mean) 
        => 1.0/2.0*(-1.0 + Math.Sqrt(1.0 + 4.0*mean*mean));

    public static double SampleParetoLim(Random rnd, double mean)
    {
        double shape = 0.5;
        double scale = GetParetoScale(mean);
        for (int i = 0; i < 1000; i++) {
            // MathNet's Pareto has support [scale, inf); shift by -scale so the
            // proportion support starts at 0, matching the shifted Pareto fit
            // (scipy loc=-xm, scale=xm) used to derive the config in lib_create_config.py.
            double sample = Pareto.Sample(rnd, scale, shape) - scale;
            if (sample <= 1)
            {
                return sample;
            }
        }
        return 1;
    }

    // For Beta(alpha, beta), mean = alpha / (alpha + beta). The configuration stores only that
    // mean; alpha stays fixed at the notebook value and beta is recovered here at runtime.
    public static double GetBetaShape(double mean, double alpha = FixedBetaAlpha)
    {
        if (mean is <= 0 or >= 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(mean), mean, "Beta mean must be strictly between zero and one.");
        }
        return alpha * (1 - mean) / mean;
    }

    public static double SampleBeta(Random rnd, double mean)
        => Beta.Sample(rnd, FixedBetaAlpha, GetBetaShape(mean));
    
    public static SexType GetSex(Random rnd, SexType sexType)
        => sexType switch
        {
            SexType.Any => rnd.CoinFlip() ? SexType.Male : SexType.Female,
            _ => sexType
        };

    public static Nucleotide SampleBase(Random rnd) 
        => (Nucleotide) rnd.Next(4);

    // Selects a single contig with probability proportional to a per-contig weight.
    // Returns null when no contig has positive weight (e.g. no contig carries a centromere).
    private static (int id, long len)? SampleContigWeighted(Random rnd, Karyotype kar, Func<int, double> weight)
    {
        var contigIds = kar.ContigIds().ToList();
        var weights = contigIds.Select(weight).ToList();
        double total = weights.Sum();
        if (total <= 0)
        {
            return null;
        }
        int idSelected = contigIds[rnd.PickRndIndex(weights, total)];
        return (idSelected, kar.ContigLen(idSelected));
    }

    private readonly record struct TerminalArm(
        int ContigId,
        long ContigLength,
        bool Direction,
        long Length);

    // Select one physical terminal arm across the whole karyotype. Sampling the candidate arms
    // directly makes both the contig and end choice proportional to usable arm length.
    private static TerminalArm? SampleTerminalArmWeighted(
        Random rnd,
        Karyotype kar,
        bool requireIntactTelomere)
    {
        var arms = kar.ContigIds()
            .SelectMany(contigId => kar.GetTerminalArms(contigId, requireIntactTelomere)
                .Select(arm => new TerminalArm(
                    contigId,
                    kar.ContigLen(contigId),
                    arm.direction,
                    arm.length)))
            .ToList();
        if (arms.Count == 0)
        {
            return null;
        }
        var weights = arms.Select(arm => (double) arm.Length).ToList();
        return arms[rnd.PickRndIndex(weights, weights.Sum())];
    }

    private static bool UsesTerminalArm(CNEventType type)
        => type is CNEventType.TailDeletion
            or CNEventType.TailDuplication
            or CNEventType.TelomereDeletion
            or CNEventType.TelomereDuplication
            or CNEventType.InternalDeletion
            or CNEventType.InternalDuplication
            or CNEventType.InternalInversion
            or CNEventType.InvertedDuplication;

    private static BaseEventData? GenerateTerminalArmEventData(
        Random rnd,
        Karyotype kar,
        CNEventPars cnEventPars)
    {
        bool requireIntactTelomere = cnEventPars.Type is
            CNEventType.TelomereDeletion or CNEventType.TelomereDuplication;
        var arm = SampleTerminalArmWeighted(rnd, kar, requireIntactTelomere);
        if (arm is not { } selected)
        {
            return null;
        }

        return cnEventPars.Type switch
        {
            CNEventType.InternalDeletion
                or CNEventType.InternalDuplication
                or CNEventType.InternalInversion
                or CNEventType.InvertedDuplication
                => new InternalEventData(
                    rnd,
                    cnEventPars,
                    selected.ContigId,
                    selected.ContigLength,
                    selected.Direction,
                    selected.Length),
            _ => new TailEventData(
                rnd,
                cnEventPars,
                selected.ContigId,
                selected.ContigLength,
                selected.Direction,
                selected.Length)
        };
    }

    // Selects contigs for events that do not use the terminal-arm path above. Other local events
    // are chosen proportional to contig length, legacy arm/centromere-bound events proportional to
    // centromere count, chromosome events uniformly among chromosome-like contigs, unrestricted
    // contig events uniformly among active contigs, and multi-contig events by random permutation.
    private static List<(int id, long len)> SelectContigs(Random rnd, Karyotype kar, CNEventType type)
    {
        switch (type)
        {
            // Multi-contig events: uniform random permutation
            case CNEventType.TIChain:
            case CNEventType.TICycle:
            case CNEventType.TIBridge:
            case CNEventType.Chromoplexy:
                return kar.ContigIds().Shuffle(rnd).Select(i => (i, kar.ContigLen(i))).ToList();
            case CNEventType.Translocation:
                return kar.ContigIds().Shuffle(rnd).Take(2).Select(i => (i, kar.ContigLen(i))).ToList();

            // Chromosome events: uniform among contigs with two intact terminal telomeres and at
            // least one centromere.
            case CNEventType.ChromDeletion:
            case CNEventType.ChromDuplication:
                return SampleContigWeighted(rnd, kar,
                    id => kar.IsChromosomeLikeContig(id) ? 1.0 : 0.0).ToList();

            // Arm / centromere-bound events: weighted by the number of centromeres in the contig
            case CNEventType.ArmDeletion:
            case CNEventType.ArmDuplication:
            case CNEventType.CentromereBoundDeletion:
            case CNEventType.CentromereBoundDuplication:
                return SampleContigWeighted(rnd, kar, id => kar.CountCentromeres(id)).ToList();

            // Within-contig events: weighted by contig length
            case CNEventType.BreakageFusionBridge:
            case CNEventType.Chromothripsis:
            case CNEventType.Pyrgo:
            case CNEventType.Rigma:
            case CNEventType.SNV:
                return SampleContigWeighted(rnd, kar, id => kar.ContigLen(id)).ToList();

            // Unrestricted contig and contig-agnostic events (Contig*, WGD, Pass, Skip): uniform pick
            default:
                return SampleContigWeighted(rnd, kar, _ => 1.0).ToList();
        }
    }

    public static BaseEventData? GenerateCNEventData(Random rnd, Karyotype kar, CNEventPars cnEventPars)
    {
        if (UsesTerminalArm(cnEventPars.Type))
        {
            return GenerateTerminalArmEventData(rnd, kar, cnEventPars);
        }

        var seq = SelectContigs(rnd, kar, cnEventPars.Type);
        if (seq.Count == 0)
        {
            return null;
        }

        switch (cnEventPars.Type)
        {
            // Whole contig events
            case CNEventType.ContigDeletion:
            case CNEventType.ContigDuplication:
            case CNEventType.ChromDeletion:
            case CNEventType.ChromDuplication:
                return new ContigEventData(cnEventPars, seq[0].id, seq[0].len);
            
            case CNEventType.Pass:
            case CNEventType.Skip:
                return new BaseEventData(cnEventPars);
            
            case CNEventType.WholeGenomeDoubling:
                return new WGDEventData(cnEventPars);
            
            // Tail-like events that are not sampled from a terminal chromosome arm.
            case CNEventType.BreakageFusionBridge:
                return new TailEventData(rnd, cnEventPars, seq[0].id, seq[0].len);

            case CNEventType.ArmDeletion:
            case CNEventType.ArmDuplication:
                return new TailEventData(rnd, cnEventPars, seq[0].id, kar.GetCentromeres(seq[0].id), seq[0].len);
            
            case CNEventType.CentromereBoundDeletion:
            case CNEventType.CentromereBoundDuplication:
                return new InternalEventData(rnd, cnEventPars, seq[0].id, kar.GetCentromeres(seq[0].id), seq[0].len);

            case CNEventType.Translocation:
                return seq.Count > 1
                    ? new PairEventData(rnd, cnEventPars, seq[0].id, seq[0].len, seq[1].id, seq[1].len)
                    : null;

            case CNEventType.Chromothripsis:
                return new ChromothripsisEventData(rnd, cnEventPars, seq[0].id, seq[0].len);

            case CNEventType.Chromoplexy:
                return new ChromoplexyEventData(rnd, cnEventPars, seq);

            case CNEventType.Pyrgo:
                return new PyrgoEventData(rnd, cnEventPars, seq[0].id, seq[0].len);

            case CNEventType.Rigma:
                return new RigmaEventData(rnd, cnEventPars, seq[0].id, seq[0].len);
            
            case CNEventType.TIChain:
            case CNEventType.TICycle:
            case CNEventType.TIBridge:
                return new TemplatedEventData(rnd, cnEventPars, seq);
            
            case CNEventType.SNV:
                return new PointMutationData(rnd, cnEventPars, seq[0].id, seq[0].len);

            default:
                throw new ArgumentOutOfRangeException(nameof(cnEventPars.Type), cnEventPars.Type, null);
        }
    }
}
