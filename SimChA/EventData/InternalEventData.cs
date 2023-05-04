using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.EventData;

public record InternalEventData : ContigEventData
{
    public long Start { get; }
    public long End { get; }

    // Constructor used for internal events
    public InternalEventData(Random rnd, Karyotype kar, CNEventP eventP, int contigId): base(eventP, contigId)
    {
        long internalSize = eventP.Get("Size", 100_000);
        long contigLen = kar.ContigLen(contigId);
        long segLen = Sampling.GetExpSeg(rnd, contigLen, internalSize);
        Start = Sampling.GetInternalPos(rnd, contigLen - segLen);
        End = Start + segLen;
    }

    public override string ToString()
    {
        return EventType switch
        {
            CNEventType.InternalDuplication => $"{EventType}\t{ContigId}\t{Start}\t{End}",
            CNEventType.InternalDeletion => $"{EventType}\t{ContigId}\t{Start}\t{End}",
            CNEventType.InternalInversion => $"{EventType}\t{ContigId}\t{Start}\t{End}",
            CNEventType.InvertedDuplication => $"{EventType}\t{ContigId}\t{Start}\t{End}",
            _ => throw new ArgumentOutOfRangeException(nameof(EventType), EventType, null)
        };
    }
    
    public override string ApplyEvent(Karyotype kar)
        => kar.ApplyEvent(this);
}
