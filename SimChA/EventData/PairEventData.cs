// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.EventData;

public record PairEventData : BaseEventData
{
    public readonly List<int>? ContigIdList = new();
    public readonly long PosA = -1;
    public readonly long PosB = -1;
    public readonly bool Direction;
    
    // Constructor used for Translocation
    public PairEventData(CNEventP eventP, List<int> contigIds, long posA, long posB, bool direction) : base(eventP)
    {
        ContigIdList = contigIds;
        PosA = posA;
        PosB = posB;
        Direction = direction;
    }
    
    public override string ApplyEvent(Karyotype kar)
        => kar.ApplyEvent(this);
}