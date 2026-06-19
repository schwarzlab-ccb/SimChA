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

    public static int SampleDiscDist(Random rnd, DistType dist, double mean)
    {
        return dist switch
        {
            DistType.Geometric => Geometric.Sample(rnd, 1 / mean),
            DistType.Poisson => Poisson.Sample(rnd, mean),
              _ => (int) mean
        };
    }

    private static double GetParetoScale(double mean) 
        => 1.0/2.0*(-1.0 + Math.Sqrt(1.0 + 4.0*mean*mean));

    public static double SampleParetoLim(Random rnd, double mean)
    {
        double shape = 0.5;
        double scale = GetParetoScale(mean);
        for (int i = 0; i < 1000; i++) {
            double sample = Pareto.Sample(rnd, scale, shape);
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

    private static List<(int id, long len)> AsList((int id, long len)? picked)
        => picked is { } p ? [p] : [];

    // Selects the contigs to be affected by the event. Within-contig events are chosen with
    // probability proportional to contig length; arm/centromere-bound events proportional to the
    // number of centromeres; multi-contig events use a uniform random permutation.
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

            // Arm / centromere-bound events: weighted by the number of centromeres in the contig
            case CNEventType.ArmDeletion:
            case CNEventType.ArmDuplication:
            case CNEventType.CentromereBoundDeletion:
            case CNEventType.CentromereBoundDuplication:
                return AsList(SampleContigWeighted(rnd, kar, id => kar.CountCentromeres(id)));

            // Within-contig events: weighted by contig length
            case CNEventType.InternalDuplication:
            case CNEventType.InternalDeletion:
            case CNEventType.InternalInversion:
            case CNEventType.InvertedDuplication:
            case CNEventType.TailDeletion:
            case CNEventType.TailDuplication:
            case CNEventType.BreakageFusionBridge:
            case CNEventType.Chromothripsis:
            case CNEventType.Pyrgo:
            case CNEventType.Rigma:
            case CNEventType.SNV:
                return AsList(SampleContigWeighted(rnd, kar, id => kar.ContigLen(id)));

            // Whole-chromosome and contig-agnostic events (Chrom*, WGD, Pass, Skip): uniform pick
            default:
                return AsList(SampleContigWeighted(rnd, kar, _ => 1.0));
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
            // Whole chromosome events
            case CNEventType.ChromDeletion:
            case CNEventType.ChromDuplication:
                return new ContigEventData(cnEventPars, seq[0].id, seq[0].len);
            
            case CNEventType.Pass:
            case CNEventType.Skip:
                return new BaseEventData(cnEventPars);
            
            case CNEventType.WholeGenomeDoubling:
                return new WGDEventData(cnEventPars);
            
            // Tail events
            case CNEventType.TailDeletion:
            case CNEventType.TailDuplication:
            case CNEventType.BreakageFusionBridge:
                return new TailEventData(rnd, cnEventPars, seq[0].id, seq[0].len);

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
