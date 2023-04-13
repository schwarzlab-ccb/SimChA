using SimChA.Simulation;

namespace SimChA.DataTypes;

public class CNEventProperties
{
    
    // Constructor used for Whole Genome Doubling
    public CNEventProperties(CNEventP eventP)
    {
        EventP = eventP;
        EventType = eventP.Type;
    }
    // Constructor used for whole-chromosome events
    public CNEventProperties(CNEventP eventP, 
        int contigId)
    {
        EventP = eventP;
        EventType = eventP.Type;
        ContigId = contigId;
    }
    // Constructor used for internal events
    public CNEventProperties(CNEventP eventP, 
        int contigId, long start, long end)
    {
        EventP = eventP;
        EventType = eventP.Type;
        ContigId = contigId;
        Start = start;
        End = end;
    }
    
    // Constructor used for Tail Events
    public CNEventProperties(CNEventP eventP,
        int contigId, long delFraction, bool delDirection)
    {
        EventP = eventP;
        EventType = eventP.Type;
        ContigId = contigId;
        DelFraction = delFraction;
        Direction = delDirection;
    }
    // Constructor used for Translocation
    public CNEventProperties(CNEventP eventP,
        List<int> contigIds, long posA, long posB, bool direction)
    {
        EventP = eventP;
        EventType = eventP.Type;
        ContigIdList = contigIds;
        PosA = posA;
        PosB = posB;
        Direction = direction;
    }

    public CNEventP EventP { get; }
    public CNEventType EventType { get; }
    public int ContigId { get; }
    public List<int>? ContigIdList { get; }
    public long Start { get; }
    public long End { get; }
    public long DelFraction { get;  }
    public bool Direction { get; }
    public long PosA { get; }
    public long PosB { get; }

    public override string ToString()
    {
        switch (EventType)
        {
            // Whole chromosome events
            case CNEventType.ChromDeletion:
            case CNEventType.ChromDuplication:
                return $"{EventType}\t{ContigId}";

            case CNEventType.WholeGenomeDoubling:
                return $"{EventType}";
            
            case CNEventType.TailDeletion:
            case CNEventType.BreakageFusionBridge:
                return $"{EventType}\t{ContigId}\t{DelFraction}\t{Direction}";
            
            case CNEventType.InternalDuplication:
            case CNEventType.InternalDeletion:
            case CNEventType.InternalInversion:
            case CNEventType.InvertedDuplication:
                return $"{EventType}\t{ContigId}\t{Start}\t{End}";
            
            case CNEventType.Translocation:
                return $"{EventType}\t{ContigIdList[0]}\t{ContigIdList[1]}\t{PosA}\t{PosB}\t{Direction}";
            default:
                throw new ArgumentOutOfRangeException(nameof(EventType), EventType, null);
        }
    }

}
