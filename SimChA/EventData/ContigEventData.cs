// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.Simulation;

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
    
    public override double GetProb()
        => 1.0;
}