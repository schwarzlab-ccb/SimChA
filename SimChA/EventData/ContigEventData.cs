using SimChA.Data;

namespace SimChA.EventData;

public record ContigEventData(CNEventPars CNEventPars, int ContigId) : BaseEventData(CNEventPars)
{
    public override void ApplyEvent(Karyotype kar)
    {
        switch (EventType)
        {
            case CNEventType.ChromDeletion:
                kar.ApplyContigDeletion(ContigId);
                break;
            case CNEventType.ChromDuplication:
                kar.ApplyContigDuplication(ContigId);
                break;
            default:
                throw new Exception($"Invalid event type {EventType} for ContigEventData");
        }
    }
    
    public override string ToString() 
        => $"contig:{ContigId}";
}