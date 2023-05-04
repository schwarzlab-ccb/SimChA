// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.EventData;

public record TemplatedEventData : BaseEventData
{
    public List<(int, long, long, bool)> Frags { get; }
    
    public TemplatedEventData(CNEventP eventP, List<(int, long, long, bool)> frags) : base(eventP)
    {
        Frags = frags;
    }
    
    public override string ToString() => $"{EventType}\t{string.Join(",", Frags)}";
    
    public override string ApplyEvent(Karyotype kar)
        => kar.ApplyEvent(this);
}