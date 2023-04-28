using SimChA.DataTypes;

namespace SimChA.EventData;

public record TailEventData : ContigEventData
{
    public readonly long DelFraction= -1;
    public readonly bool Direction;
    
    // Constructor used for Tail Events
    public TailEventData(CNEventP eventP, int contigId, long delFraction, bool delDirection) : base(eventP, contigId)
    {
        DelFraction = delFraction;
        Direction = delDirection;
    }

    public override string ToString()
    {
        return EventType switch
        {
            CNEventType.TailDeletion => $"{EventType}\t{ContigId}\t{DelFraction}\t{Direction}",
            CNEventType.BreakageFusionBridge => $"{EventType}\t{ContigId}\t{DelFraction}\t{Direction}",
            _ => throw new ArgumentOutOfRangeException(nameof(EventType), EventType, null)
        };
    }

}
