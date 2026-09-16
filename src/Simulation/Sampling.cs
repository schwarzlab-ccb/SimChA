using System.Collections.Concurrent;
using SimChA.Computation;
using SimChA.EventData;
using MathNet.Numerics.Distributions;
using SimChA.Data;

namespace SimChA.Simulation;

public static class Sampling
{
    public const double FixedBetaAlpha = 0.6;
    public const double FixedParetoShape = 0.5;

    // The bounded Pareto flattens toward uniform on [0,1] as its scale grows, whatever the shape,
    // so no shape can reach a mean of one half. Scales are searched over this bracket, which spans
    // means from far below anything a config carries up to 0.49999.
    private const double ParetoScaleMin = 1e-300;
    private const double ParetoScaleMax = 1e4;

    // Solving for the scale costs a bisection, and the (mean, shape) pair is fixed per event type
    // for a whole run, so it is solved once rather than per event.
    private static readonly ConcurrentDictionary<(double Mean, double Shape), double> ParetoScaleCache = new();

    public static long GetNormSeg(Random rnd, long contigLen, double meanFrac) 
        => (long) Math.Clamp(Math.Round(contigLen * Normal.Sample(rnd, meanFrac, meanFrac / 3)), 1, contigLen);
    
    public static long GetExpSeg(Random rnd, long contigLen, long meanLen) 
        => (long) Math.Clamp(Math.Round(meanLen * Exponential.Sample(rnd, 1)), 1, contigLen);
    
    public static long GetExpSeg(Random rnd, long contigLen, double meanFrac)
        => (long) Math.Clamp(Math.Round(contigLen * meanFrac * Exponential.Sample(rnd, 1)), 1, contigLen);

    public static long GetParetoSeg(Random rnd, long armLength, double meanFrac, double? shape = null)
    {
        // The same saturating guards GetBetaSeg applies. A mean at or above the whole arm asks for
        // the longest event available, which is what the rejection loop this replaced effectively
        // returned: at a scale that large every one of its thousand tries exceeded the arm and it
        // fell through to a full-length draw. A mean between one half and one is a different case
        // and is rejected below, because no bounded Pareto has a mean there.
        if (meanFrac >= 1)
        {
            return armLength;
        }
        if (meanFrac <= 0)
        {
            return 1;
        }
        return (long) Math.Clamp(
            Math.Round(armLength * SampleParetoLim(rnd, meanFrac, shape)),
            1,
            armLength);
    }

    public static long GetBetaSeg(Random rnd, long armLength, double meanFrac, double? alpha = null)
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
            Math.Round(armLength * SampleBeta(rnd, meanFrac, alpha)),
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

    // log(1 + x) and exp(x) - 1, accurate for a small argument. .NET's double.LogP1 and
    // double.ExpM1 are the naive Math.Log(1 + x) / Math.Exp(x) - 1 and lose most of their
    // significance below about 1e-9 -- at 1e-12 they are wrong in the fifth digit, which is enough
    // to make BoundedParetoMean return 1.9e9 for a quantity bounded by one half.
    private static double Log1P(double x)
    {
        double sum = 1.0 + x;
        return sum == 1.0 ? x : Math.Log(sum) * x / (sum - 1.0);
    }

    private static double ExpM1(double x)
    {
        double exp = Math.Exp(x);
        if (exp == 1.0)
        {
            return x;
        }
        return exp - 1.0 == -1.0 ? -1.0 : (exp - 1.0) * x / Math.Log(exp);
    }

    // Mean of the bounded shifted Pareto on [0,1]: E[X - L | X <= L + 1] for X ~ Pareto(L, a).
    public static double BoundedParetoMean(double scale, double shape)
    {
        double logRatio = Log1P(1.0 / scale);
        // 1 - (L/(L+1))^a, formed directly rather than by subtracting a near-one quantity from one.
        double below = -ExpM1(-shape * logRatio);
        double integral;
        if (Math.Abs(shape - 1) < 1e-9)
        {
            integral = scale * logRatio;
        }
        else if (scale > 1)
        {
            // Above one, (L+1)^(1-a) and L^(1-a) agree to many digits, so the direct difference
            // cancels to noise -- it returns means above 300 at L = 1e9. The identity
            // L^a*((L+1)^(1-a) - L^(1-a)) == L*expm1((1-a)*log((L+1)/L)) avoids the subtraction.
            integral = scale * ExpM1((1 - shape) * logRatio) / (1 - shape);
        }
        else
        {
            // Below one the two powers are far apart and the direct form is exact, while the
            // expm1 form would overflow: log((L+1)/L) grows without bound as the scale shrinks,
            // reaching 691 at the smallest scale searched. L^a is distributed into the bracket
            // because L^a*L^(1-a) is just L, whereas evaluating the two powers separately gives
            // 0 * infinity at a scale of 1e-300 once the shape passes one.
            integral = (Math.Pow(scale, shape) * Math.Pow(scale + 1, 1 - shape) - scale)
                       / (1 - shape);
        }
        return (integral - (1 - below)) / below;
    }

    // The scale that makes the bounded distribution's mean equal the configured Frac, at this shape.
    // The formula this replaced, 1/2*(-1 + sqrt(1 + 4*mean^2)), contained no shape at all: it was
    // not the inverse of any mean, and left the realized proportion 10-13% below the configured one.
    public static double GetParetoScale(double mean, double shape = FixedParetoShape)
    {
        if (shape <= 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(shape), shape, "Pareto shape must be positive.");
        }
        double lowest = BoundedParetoMean(ParetoScaleMin, shape);
        double highest = BoundedParetoMean(ParetoScaleMax, shape);
        if (mean <= lowest || mean >= highest)
        {
            throw new ArgumentOutOfRangeException(
                nameof(mean), mean,
                $"At shape {shape:G6} the bounded Pareto mean must lie in "
                + $"({lowest:G6}, {highest:G6}); it tends to uniform on [0,1] and so cannot reach "
                + "one half.");
        }
        // BoundedParetoMean is strictly increasing in the scale, so plain bisection converges; in
        // log space 200 halvings take the bracket far below double precision.
        double lo = Math.Log(ParetoScaleMin), hi = Math.Log(ParetoScaleMax);
        for (int i = 0; i < 200; i++)
        {
            double mid = 0.5 * (lo + hi);
            if (BoundedParetoMean(Math.Exp(mid), shape) < mean)
            {
                lo = mid;
            }
            else
            {
                hi = mid;
            }
        }
        return Math.Exp(0.5 * (lo + hi));
    }

    // Bounded shifted Pareto on [0,1], drawn by inverse CDF. The previous version sampled the
    // unbounded distribution and rejected anything above 1, which is the same distribution but
    // costs 1.2 draws per sample at shape 0.5 and 2.0 at shape 0.12, and needed a fallback for the
    // case where every try was rejected. The quantile is exact, so one draw always suffices.
    public static double SampleParetoLim(Random rnd, double mean, double? shape = null)
    {
        double paretoShape = shape ?? FixedParetoShape;
        double scale = ParetoScaleCache.GetOrAdd(
            (mean, paretoShape), key => GetParetoScale(key.Mean, key.Shape));
        double below = -ExpM1(-paretoShape * Log1P(1.0 / scale));   // 1 - (L/(L+1))^a
        double sample = scale * Math.Pow(1 - rnd.NextDouble() * below, -1.0 / paretoShape) - scale;
        // NextDouble draws from [0,1), so the quantile lands in [0,1); the clamp only guards the
        // rounding at the top end, where 1 - u*below can underflow at an extreme shape.
        return Math.Clamp(sample, 0, 1);
    }

    // For Beta(alpha, beta), mean = alpha / (alpha + beta). The configuration stores the mean and
    // optionally the alpha; beta is recovered here at runtime so the mean comes out exactly.
    public static double GetBetaShape(double mean, double alpha = FixedBetaAlpha)
    {
        if (alpha <= 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(alpha), alpha, "Beta alpha must be positive.");
        }
        if (mean is <= 0 or >= 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(mean), mean, "Beta mean must be strictly between zero and one.");
        }
        return alpha * (1 - mean) / mean;
    }

    public static double SampleBeta(Random rnd, double mean, double? alpha = null)
    {
        double betaAlpha = alpha ?? FixedBetaAlpha;
        return Beta.Sample(rnd, betaAlpha, GetBetaShape(mean, betaAlpha));
    }
    
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

    private readonly record struct ContigArm(int ContigId, long ContigLength, TerminalArm Arm);

    // Select one physical terminal arm across the whole karyotype. Sampling the candidate arms
    // directly makes both the contig and end choice proportional to usable arm length.
    private static ContigArm? SampleTerminalArmWeighted(
        Random rnd,
        Karyotype kar,
        bool requireIntactTelomere)
    {
        var arms = kar.ContigIds()
            .SelectMany(contigId => kar.GetTerminalArms(contigId, requireIntactTelomere)
                .Select(arm => new ContigArm(contigId, kar.ContigLen(contigId), arm)))
            .ToList();
        if (arms.Count == 0)
        {
            return null;
        }
        var weights = arms.Select(arm => (double) arm.Arm.UsableLength).ToList();
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
                    rnd, cnEventPars, selected.ContigId, selected.ContigLength, selected.Arm),
            _ => new TailEventData(
                rnd, cnEventPars, selected.ContigId, selected.ContigLength, selected.Arm)
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
                return kar.ContigIds()
                    .Where(id => kar.CountCentromeres(id) == 1 &&
                                 kar.GetTerminalArms(id, false).Count > 0)
                    .Shuffle(rnd).Take(2).Select(i => (i, kar.ContigLen(i))).ToList();

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

            // BFB requires a monocentric contig with an arm available for the initial break.
            case CNEventType.BreakageFusionBridge:
                return SampleContigWeighted(rnd, kar,
                    id => kar.CountCentromeres(id) == 1 &&
                          kar.GetTerminalArms(id, false).Count > 0
                        ? kar.ContigLen(id)
                        : 0).ToList();

            // Other within-contig events: weighted by contig length
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
            
            case CNEventType.BreakageFusionBridge:
            {
                var arms = kar.GetTerminalArms(seq[0].id, false);
                var arm = arms[rnd.Next(arms.Count)];
                return new BFBEventData(rnd, cnEventPars, seq[0].id, seq[0].len, arm);
            }

            case CNEventType.ArmDeletion:
            case CNEventType.ArmDuplication:
                return new TailEventData(rnd, cnEventPars, seq[0].id, kar.GetCentromeres(seq[0].id), seq[0].len);
            
            case CNEventType.CentromereBoundDeletion:
            case CNEventType.CentromereBoundDuplication:
                return new InternalEventData(rnd, cnEventPars, seq[0].id, kar.GetCentromeres(seq[0].id), seq[0].len);

            case CNEventType.Translocation:
            {
                if (seq.Count < 2)
                {
                    return null;
                }
                var armsA = kar.GetTerminalArms(seq[0].id, false);
                var armsB = kar.GetTerminalArms(seq[1].id, false);
                var armA = armsA[rnd.Next(armsA.Count)];
                var armB = armsB[rnd.Next(armsB.Count)];
                return new PairEventData(
                    rnd, cnEventPars,
                    seq[0].id, seq[0].len, armA,
                    seq[1].id, seq[1].len, armB);
            }

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
