// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.EventData;

public record BaseEventData(CNEventPars CNEventPars)
{
    public CNEventType EventType => CNEventPars.Type;
    
    public virtual void ApplyEvent(Karyotype kar)
    {
        kar.ApplyWGD();
    }

    public override string ToString()
        => "";
}