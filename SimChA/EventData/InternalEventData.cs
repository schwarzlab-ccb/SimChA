using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.EventData;

public record InternalEventData : ContigEventData
{
    public readonly long Start = -1;
    public readonly long End = -1;

    // Constructor used for internal events
    public InternalEventData(CNEventP eventP, int contigId, long start, long end) : base(eventP, contigId)
    {
        Start = start;
        End = end;
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
