using System.ComponentModel.DataAnnotations;
using Extreme.Statistics.Distributions;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;

namespace SimChA.Simulation;

public static class Sampling
{
    public static long GetNormSeg(Random rnd, long contigLen, double meanFrac) 
        => Math.Max(1, Math.Min((long) Math.Round(contigLen * NormalDistribution.Sample(rnd, meanFrac, meanFrac / 3)), contigLen));
    
    public static long GetExpSeg(Random rnd, long contigLen, long meanLen) 
        => Math.Max(1, Math.Min((long) Math.Round(ExponentialDistribution.Sample(rnd, meanLen)), contigLen));
    
    public static long GetExpSeg(Random rnd, long contigLen, double meanFrac) 
        => Math.Max(1, Math.Min((long) Math.Round(contigLen * ExponentialDistribution.Sample(rnd, meanFrac)), contigLen));
    
    public static double GetExpProb(long segLen, double scale)
        => Math.Exp(-segLen / scale) / scale;
    
    public static long GetPos(Random rnd, long contigLen)
        => rnd.NextInt64(0, contigLen);
    
    public static int GetFragCount(Random rnd, double mean)
        => (int) Math.Max(1, Math.Round(NormalDistribution.Sample(rnd, mean)));

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
    
    public static List<double> CreateRandomMixture(Random rnd, double[] concentrations)
        => concentrations.Any() ? new DirichletDistribution(concentrations).Sample(rnd).ToList() : new List<double>();
    
    public static double SampleDist(Random rnd, DataTypes.Distribution dist, double mean = 1)
    {
        return dist switch
        {
            DataTypes.Distribution.Exponential => ExponentialDistribution.Sample(rnd, mean),
            DataTypes.Distribution.Normal => NormalDistribution.Sample(rnd, mean, .5),
            _ => 1
        };
    }

    public static int SampleDistance(Random rnd, DataTypes.Distribution dist, double mean)
    {
        return dist switch
        {
            DataTypes.Distribution.Exponential => (int) Math.Round( ExponentialDistribution.Sample(rnd, mean)),
            DataTypes.Distribution.Geometric => GeometricDistribution.Sample(rnd, mean),
            DataTypes.Distribution.Poisson => PoissonDistribution.Sample(rnd, mean),
            DataTypes.Distribution.Normal => throw new Exception("Normal distribution not supported for distance sampling"),
            _ => (int) mean
        };
    }

    public static int SampleDistInt(Random rnd, DataTypes.Distribution dist, double mean)
    {
        return dist switch
        {
            DataTypes.Distribution.Geometric => GeometricDistribution.Sample(rnd, mean),
            DataTypes.Distribution.Poisson => PoissonDistribution.Sample(rnd, mean),
            DataTypes.Distribution.Unit => (int)mean,
            _ =>1
        };
    }

    public static SexEnum GetSex(Random rnd, SexEnum sexEnum)
        => sexEnum switch
        {
            SexEnum.None => rnd.CoinFlip() ? SexEnum.Male : SexEnum.Female,
            _ => sexEnum
        };

    public static Nucleotide SampleBase(Random rnd) 
        => (Nucleotide) rnd.Next(4);

    private static (int id, long len) SampleContigsByLength(Random rnd, Karyotype kar)
    {
        // Karyotype stores 0-length contigs for contig-ID-preservation, so we need to filter them out
        var contigIds = kar.ContigIds().Where(i => kar.ContigLen(i) > 0).ToList();
        long totalLength = contigIds.Sum(kar.ContigLen);
        var pArray = contigIds.Select(i => kar.ContigLen(i)/(1.0*totalLength)).ToList();
        int idSelected = contigIds.ToList()[rnd.PickRndIndex(pArray)];
        return (idSelected, kar.ContigLen(idSelected));
    }
    
    public static BaseEventData? GenerateCNEventData(Random rnd, Karyotype kar, CNEventPars cnEventPars)
    {
        List<(int id, long len)> seq = kar.ContigIds().Shuffle(rnd).Select(i => (i, kar.ContigLen(i))).ToList();
        if (seq.Count == 0)
            return null;
        
        switch (cnEventPars.Type)
        {
            // Whole chromosome events
            case CNEventType.ChromDeletion:
            case CNEventType.ChromDuplication:
                return new ContigEventData(cnEventPars, seq[0].id);
            
            case CNEventType.Pass:
            case CNEventType.WholeGenomeDoubling:
                return new BaseEventData(cnEventPars);
            
            // Tail events
            case CNEventType.TailDeletion:
            case CNEventType.TailDuplication:
            case CNEventType.BreakageFusionBridge:
                return new TailEventData(rnd, cnEventPars, seq[0].id, seq[0].len);

            case CNEventType.ArmDeletion:
            case CNEventType.ArmDuplication:
                var cents = kar.GetCentromeres(seq[0].id);
                return cents.Count == 0 
                    ? null
                    : new TailEventData(rnd, cnEventPars, seq[0].id, cents);

            case CNEventType.CentromereBoundDeletion:
            case CNEventType.CentromereBoundDuplication:
                cents = kar.GetCentromeres(seq[0].id);
                return cents.Count == 0
                    ? null
                    : new InternalEventData(rnd, cnEventPars, seq[0].id, seq[0].len, cents);

            // Internal events
            case CNEventType.InternalDuplication:
            case CNEventType.InternalDeletion:
            case CNEventType.InternalInversion:
            case CNEventType.InvertedDuplication:
                return new InternalEventData(rnd, cnEventPars, seq[0].id, seq[0].len);
            
            case CNEventType.Translocation:
                return seq.Count < 2 
                    ? null 
                    : new PairEventData(rnd, cnEventPars, seq[0].id, seq[0].len, seq[1].id, seq[0].len);

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
                var pointSeq = SampleContigsByLength(rnd, kar);
                return new PointMutationData(rnd, cnEventPars, pointSeq.id, pointSeq.len);

            default:
                throw new ArgumentOutOfRangeException(nameof(cnEventPars.Type), cnEventPars.Type, null);
        }
    }
}
