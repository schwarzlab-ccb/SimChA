// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using Extreme.Mathematics;
using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.EventData;

public record RigmaEventData : BaseEventData
{
    public readonly int ContigId = -1;
    public readonly long Start = -1;
    public readonly List<long>? StopsList = new();
    public RigmaEventData(CNEventP eventP, int contigId, long startPoint, List<long> stops) : base(eventP)
    {
        ContigId = contigId;
        Start = startPoint;
        StopsList = stops;
    }
    
    public override string ApplyEvent(Karyotype kar)
        => kar.ApplyEvent(this);
}
