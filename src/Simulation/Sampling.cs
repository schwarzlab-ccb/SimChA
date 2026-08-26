using SimChA.Computation;
using SimChA.EventData;
using MathNet.Numerics.Distributions;
using SimChA.Data;

namespace SimChA.Simulation;

public static class Sampling
{
    public static long GetNormSeg(Random rnd, long contigLen, double meanFrac) 
        => (long) Math.Clamp(Math.Round(contigLen * Normal.Sample(rnd, meanFrac, meanFrac / 3)), 1, contigLen);
    
    public static long GetExpSeg(Random rnd, long contigLen, long meanLen) 
        => (long) Math.Clamp(Math.Round(meanLen * Exponential.Sample(rnd, 1)), 1, contigLen);
    
    public static long GetExpSeg(Random rnd, long contigLen, double meanFrac)
        => (long) Math.Clamp(Math.Round(contigLen * meanFrac * Exponential.Sample(rnd, 1)), 1, contigLen);

    public static long GetParetoSeg(Random rnd, long contigLen, double meanFrac)
        => (long) Math.Clamp(Math.Round(contigLen * SampleParetoLim(rnd, meanFrac)), 1, contigLen);

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

    // Selects the contigs to be affected by the event. Within-contig events are chosen with
    // probability proportional to contig length; arm/centromere-bound events proportional to the
    // number of centromeres; telomere events proportional to eligible terminal telomeres;
    // chromosome events uniformly among chromosome-like contigs; unrestricted contig events
    // uniformly among active contigs; and multi-contig events use a uniform random permutation.
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

            // Tail loss retains length weighting but excludes contigs without an inward centromere.
            case CNEventType.TailDeletion:
                return SampleContigWeighted(rnd, kar,
                    id => kar.CountTailLossEnds(id) > 0 ? kar.ContigLen(id) : 0).ToList();

            // Telomere loss: weighted by intact terminal telomeres that have a centromere on their
            // inward side. The distance to the nearest such centromere bounds the deletion.
            case CNEventType.TelomereDeletion:
                return SampleContigWeighted(rnd, kar, id => kar.CountTelomereBoundLossEnds(id)).ToList();

            // Telomere gain: weighted by the number of intact telomeres at physical contig ends.
            case CNEventType.TelomereDuplication:
                return SampleContigWeighted(rnd, kar, id => kar.CountIntactTelomereEnds(id)).ToList();

            // Within-contig events: weighted by contig length
            case CNEventType.InternalDuplication:
            case CNEventType.InternalDeletion:
            case CNEventType.InternalInversion:
            case CNEventType.InvertedDuplication:
            case CNEventType.TailDuplication:
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
            
            // Tail events
            case CNEventType.TailDuplication:
            case CNEventType.BreakageFusionBridge:
                return new TailEventData(rnd, cnEventPars, seq[0].id, seq[0].len);

            case CNEventType.TailDeletion:
                var tailLossDirections = kar.GetTailLossDirections(seq[0].id);
                if (tailLossDirections.Count == 0)
                {
                    return null;
                }
                var tailLossDirection = tailLossDirections[rnd.Next(tailLossDirections.Count)];
                return new TailEventData(
                    rnd,
                    cnEventPars,
                    seq[0].id,
                    seq[0].len,
                    tailLossDirection.direction,
                    tailLossDirection.maxLength);

            case CNEventType.TelomereDuplication:
                var directions = kar.GetIntactTelomereDirections(seq[0].id);
                return directions.Count == 0
                    ? null
                    : new TailEventData(
                        rnd,
                        cnEventPars,
                        seq[0].id,
                        seq[0].len,
                        directions[rnd.Next(directions.Count)]);

            case CNEventType.TelomereDeletion:
                var lossDirections = kar.GetTelomereBoundLossDirections(seq[0].id);
                if (lossDirections.Count == 0)
                {
                    return null;
                }
                var lossDirection = lossDirections[rnd.Next(lossDirections.Count)];
                return new TailEventData(
                    rnd,
                    cnEventPars,
                    seq[0].id,
                    seq[0].len,
                    lossDirection.direction,
                    lossDirection.maxLength);

            case CNEventType.ArmDeletion:
            case CNEventType.ArmDuplication:
                return new TailEventData(rnd, cnEventPars, seq[0].id, kar.GetCentromeres(seq[0].id), seq[0].len);
            
            case CNEventType.CentromereBoundDeletion:
            case CNEventType.CentromereBoundDuplication:
                return new InternalEventData(rnd, cnEventPars, seq[0].id, kar.GetCentromeres(seq[0].id), seq[0].len);

            // Internal events
            case CNEventType.InternalDuplication:
            case CNEventType.InternalDeletion:
            case CNEventType.InternalInversion:
            case CNEventType.InvertedDuplication:
                return new InternalEventData(rnd, cnEventPars, seq[0].id, seq[0].len);
            
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
