// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using Extreme.Mathematics;
using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.EventData;

public record PyrgoEventData : BaseEventData
{
    public readonly int ContigId = -1;
    public readonly List<(long, long)>? FragmentsList = new();
    public PyrgoEventData(CNEventP eventP, int contigId, List<(long,long)> frags) : base(eventP)
    {
        ContigId = contigId;
        FragmentsList = frags;
    }
    public override string ToString() => $"{EventType}\t{ContigId}";
    
    public override string ApplyEvent(Karyotype kar)
        => kar.ApplyEvent(this);
}