using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.EventData;
using MathNet.Numerics.Distributions;
using SimChA.Data;

namespace SimChA.Simulation;

public static class Sampling
{
    public static long GetNormSeg(Random rnd, long contigLen, double meanFrac) 
        => Math.Max(1, Math.Min((long) Math.Round(contigLen * Normal.Sample(rnd, meanFrac, meanFrac / 3)), contigLen));
    
    public static long GetExpSeg(Random rnd, long contigLen, long meanLen) 
        => Math.Max(1, Math.Min((long) Math.Round(Exponential.Sample(rnd, meanLen)), contigLen));
    
    public static long GetExpSeg(Random rnd, long contigLen, double meanFrac) 
        => Math.Max(1, Math.Min((long) Math.Round(contigLen * Exponential.Sample(rnd, meanFrac)), contigLen));
    
    public static double GetExpProb(long segLen, double scale)
        => Math.Exp(-segLen / scale) / scale;
    
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
    
    public static List<double> CreateRandomMixture(Random rnd, double[] concentrations)
        => concentrations.Any() ? new Dirichlet(concentrations, rnd).Sample().ToList() : new List<double>();
    
    public static double SampleDist(Random rnd, DistType dist, double mean = 1)
    {
        return dist switch
        {
            DistType.Exponential => Exponential.Sample(rnd, 1 / mean),
            DistType.Normal => Normal.Sample(rnd, mean, .5),
            DistType.Geometric => Geometric.Sample(rnd, 1 / mean),
            DistType.Poisson => Poisson.Sample(rnd, mean),
            _ => mean
        };
    }

    public static int SampleDistInt(Random rnd, DistType dist, double mean)
    {
        return dist switch
        {
            DistType.Geometric => Geometric.Sample(rnd, 1 / mean),
            DistType.Poisson => Poisson.Sample(rnd, mean),
            DistType.Normal => throw new Exception($"{dist} distribution not supported for distance sampling"),
            DistType.Exponential => throw new Exception($"{dist} distribution not supported for distance sampling"),
            _ => (int) mean
        };
    }
    
    public static SexType GetSex(Random rnd, SexType sexType)
        => sexType switch
        {
            SexType.None => rnd.CoinFlip() ? SexType.Male : SexType.Female,
            _ => sexType
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
