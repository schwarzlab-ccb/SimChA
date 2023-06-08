using SimChA.Simulation;

namespace SimChA.EventData;

public record InternalEventData : ContigEventData
{
    public long Start { get; }
    public long End { get; }

    // Constructor used for internal events
    public InternalEventData(Random rnd, CNEventPars CNEventPars, int contigId, long contigLen) : base(CNEventPars, contigId)
    {
        long internalSize = CNEventPars.GetInt("Size");
        long segLen = Sampling.GetExpSeg(rnd, contigLen, internalSize);
        Start = Sampling.GetInternalPos(rnd, contigLen - segLen);
        End = Start + segLen;
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
            kar.ApplyInternalDuplication(ContigId, Start, End);
        }
    }

    public override string ToString()
        => $"contig:{ContigId};start:{Start};end:{End}";
}