using System.Security.Cryptography;
using SimChA.Computation;
using SimChA.Simulation;

namespace SimChA.EventData;

public record InternalEventData : ContigEventData
{
    public long Start { get; }
    public long End { get; }

    // Constructor used for internal events
    public InternalEventData(Random rnd, CNEventPars CNEventPars, int contigId, long contigLen) : base(CNEventPars, contigId)
    {
        long segLen = Sampling.GetExpSeg(rnd, contigLen, CNEventPars.Size);
        //Start = rnd.NextInt64(segLen, contigLen-segLen);
        Start = Sampling.GetPos(rnd, contigLen - segLen); 
        End = Start + segLen;
    }

    // Constructor for Centromere-bound events
    public InternalEventData(Random rnd, CNEventPars CNEventPars, int contigId, long contigLen,
        IEnumerable<(long start, long end)> centromeres) : base(CNEventPars, contigId)
    {
        var cent = centromeres.Shuffle(rnd).First();
        long segLen = Sampling.GetExpSeg(rnd, contigLen, CNEventPars.Size);
        var pos = rnd.NextInt64(cent.start, cent.end);
        if (rnd.CoinFlip())
        {
            Start = Math.Max(0, pos);
            End = Math.Min(pos + segLen, contigLen);
        }
        else
        {
            Start = Math.Max(0, pos - segLen);
            End = Math.Min(pos, contigLen);
        }
    }

    public override void ApplyEvent(Karyotype kar)
    {
        switch (EventType)
        {
            case CNEventType.InternalDuplication:
            case CNEventType.CentromereBoundDuplication:
                kar.ApplyInternalDuplication(ContigId, Start, End);
                break;
            case CNEventType.InternalDeletion:
            case CNEventType.CentromereBoundDeletion:
                kar.ApplyInternalDeletion(ContigId, Start, End);
                break;
            case CNEventType.InternalInversion:
                kar.ApplyInternalInversion(ContigId, Start, End);
                break;
            case CNEventType.InvertedDuplication:
                kar.ApplyInvertedDuplication(ContigId, Start, End);
                break;
            default:
                throw new Exception($"Invalid event type {EventType} for InternalEventData");
        }
    }

    public override string ToString()
        => $"contig:{ContigId};start:{Start};end:{End}";
}