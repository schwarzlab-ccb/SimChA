// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.EventData;

public record BaseEventData(CNEventP EventP)
{
    public CNEventType EventType => EventP.Type;
    public override string ToString() => $"{EventType}";
    public virtual string ApplyEvent(Karyotype kar)
        => kar.ApplyEvent(this);
}