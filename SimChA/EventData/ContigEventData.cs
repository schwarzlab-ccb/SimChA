// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.EventData;

public record ContigEventData(CNEventPars CNEventPars, int ContigId) : BaseEventData(CNEventPars)
{
    public override void ApplyEvent(Karyotype kar)
    {
        if (EventType == CNEventType.ChromDeletion)
        {
            kar.ApplyContigDeletion(ContigId);
        }
        else if (EventType == CNEventType.ChromDuplication)
        {
            kar.ApplyContigDuplication(ContigId);
        }
    }
    
    public override string ToString() 
        => $"contig:{ContigId}";
}