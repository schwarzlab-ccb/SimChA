using System.Security.Cryptography;
using SimChA.Simulation;

namespace SimChA.EventData;

public record InternalEventData : ContigEventData
{
    public long Start { get; }
    public long End { get; }
    public double Prob { get; }

    // Constructor used for internal events
    public InternalEventData(Random rnd, CNEventPars CNEventPars, int contigId, long contigLen) : base(CNEventPars, contigId)
    {
        long segLen = Sampling.GetExpSeg(rnd, contigLen, CNEventPars.Size);
        Start = Sampling.GetPos(rnd, contigLen - segLen);
        End = Start + segLen;
        Prob = Sampling.GetExpProb(segLen, CNEventPars.Size);
    }

    public override void ApplyEvent(Karyotype kar)
    {
        if (EventType == CNEventType.InternalDuplication)
        {
            kar.ApplyInternalDuplication(ContigId, Start, End);
        }
        else if (EventType == CNEventType.InternalDeletion)
        {
            kar.ApplyInternalDeletion(ContigId, Start, End);
        }
        else if (EventType == CNEventType.InternalInversion)
        {
            kar.ApplyInternalInversion(ContigId, Start, End);
        }
        else if (EventType == CNEventType.InvertedDuplication)
        {
            kar.ApplyInvertedDuplication(ContigId, Start, End);
        }
    }

    public override string ToString()
        => $"contig:{ContigId};start:{Start};end:{End}";
    public override double GetProb()
        => Prob;
}