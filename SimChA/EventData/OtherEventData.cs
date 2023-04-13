using SimChA.DataTypes;

namespace SimChA.EventData;

public record OtherEventData : ContigEventData
{
    public readonly long Start = -1;
    public readonly long End = -1;
    public readonly long DelFraction= -1;
    public readonly bool Direction;

    // Constructor used for internal events
    public OtherEventData(CNEventP eventP, int contigId, long start, long end) : base(eventP, contigId)
    {
        Start = start;
        End = end;
    }
    
    // Constructor used for Tail Events
    public OtherEventData(CNEventP eventP, int contigId, long delFraction, bool delDirection) : base(eventP, contigId)
    {
        DelFraction = delFraction;
        Direction = delDirection;
    }

    public override string ToString()
    {
        return EventType switch
        {
            // Whole chromosome events
            CNEventType.ChromDeletion => $"{EventType}\t{ContigId}",
            CNEventType.ChromDuplication => $"{EventType}\t{ContigId}",
            CNEventType.WholeGenomeDoubling => $"{EventType}",
            CNEventType.TailDeletion => $"{EventType}\t{ContigId}\t{DelFraction}\t{Direction}",
            CNEventType.BreakageFusionBridge => $"{EventType}\t{ContigId}\t{DelFraction}\t{Direction}",
            CNEventType.InternalDuplication => $"{EventType}\t{ContigId}\t{Start}\t{End}",
            CNEventType.InternalDeletion => $"{EventType}\t{ContigId}\t{Start}\t{End}",
            CNEventType.InternalInversion => $"{EventType}\t{ContigId}\t{Start}\t{End}",
            CNEventType.InvertedDuplication => $"{EventType}\t{ContigId}\t{Start}\t{End}",
            _ => throw new ArgumentOutOfRangeException(nameof(EventType), EventType, null)
        };
    }

}
